"""Agent-facing tool façade for pyHMT2D workflows."""

import csv
import hashlib
import json
import os
import random
import shutil
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, List

import pyHMT2D

from .errors import JobExecutionError, ValidationError
from .models import ArtifactRecord, ToolResponse, utc_now_iso
from .schemas import (
    validate_calibration_request,
    validate_convert_request,
    validate_modify_srh_case_request,
    validate_monte_carlo_srh_request,
    validate_result_summary_request,
    validate_run_model_request,
)


@contextmanager
def _pushd(path: Path):
    original = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(original)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _artifact(path: Path, with_hash: bool = False) -> ArtifactRecord:
    exists = path.exists()
    size_bytes = path.stat().st_size if exists else None
    sha = _sha256(path) if exists and with_hash and path.is_file() else None
    return ArtifactRecord(path=str(path), exists=exists, size_bytes=size_bytes, sha256=sha)


def write_run_manifest(output_dir: Path, request: Dict[str, Any], response: ToolResponse) -> Path:
    """Write a run manifest for observability and reproducibility."""
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "run_manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(
            {
                "created_at": utc_now_iso(),
                "request": request,
                "response": response.to_dict(),
                "cwd": os.getcwd(),
            },
            f,
            indent=2,
        )
    return manifest_path


def _save_modified_srh_control(data_obj: Any, output_control_file: str | None) -> Path:
    if data_obj.control_type == "SRHHydro":
        if output_control_file:
            data_obj.srhhydro_obj.save_as(output_control_file)
            return Path(output_control_file)

        data_obj.srhhydro_obj.save_as()
        return Path(data_obj.srhhydro_obj.srhhydro_filename)

    if output_control_file:
        data_obj.srhsif_obj.save_as(output_control_file)
        return Path(output_control_file)

    data_obj.srhsif_obj.save_as()
    return Path(data_obj.srhsif_obj.srhsif_filename)


def validate_job_spec(payload: Dict[str, Any]) -> ToolResponse:
    """Validate a typed job spec without executing solvers."""
    response = ToolResponse(status="succeeded", message="Validation successful")

    try:
        request_type = payload.get("request_type")
        if request_type == "run_model":
            normalized = validate_run_model_request(payload)
        elif request_type == "calibration":
            normalized = validate_calibration_request(payload)
        elif request_type == "convert":
            normalized = validate_convert_request(payload)
        elif request_type == "result_summary":
            normalized = validate_result_summary_request(payload)
        elif request_type == "modify_srh_case":
            normalized = validate_modify_srh_case_request(payload)
        elif request_type == "monte_carlo_srh":
            normalized = validate_monte_carlo_srh_request(payload)
        else:
            raise ValidationError(
                "request_type must be one of: run_model, calibration, convert, result_summary, "
                "modify_srh_case, monte_carlo_srh"
            )

        response.metrics["normalized_request"] = normalized
    except ValidationError as exc:
        response.status = "failed"
        response.error_code = "VALIDATION_ERROR"
        response.message = str(exc)

    return response.finalize()


def run_calibration(payload: Dict[str, Any]) -> ToolResponse:
    """Run calibration job from a validated configuration JSON path."""
    response = ToolResponse(status="running", message="Calibration started")

    try:
        req = validate_calibration_request(payload)
        config_path = Path(req["config_json"])

        response.metrics["config_json"] = str(config_path)
        if req["dry_run"]:
            response.status = "succeeded"
            response.message = "Calibration dry-run validation complete"
        else:
            calibrator = pyHMT2D.Calibration.Calibrator(str(config_path))
            calibrator.calibrate()
            response.status = "succeeded"
            response.message = "Calibration completed"
            response.artifacts.append(_artifact(Path("calibration.log"), with_hash=True))

        manifest = write_run_manifest(config_path.parent, req, response)
        response.artifacts.append(_artifact(manifest, with_hash=True))
    except ValidationError as exc:
        response.status = "failed"
        response.error_code = "VALIDATION_ERROR"
        response.message = str(exc)
    except SystemExit as exc:
        response.status = "failed"
        response.error_code = "MODEL_EXIT"
        response.message = f"Model exited during calibration: {exc}"
    except Exception as exc:  # noqa: BLE001
        response.status = "failed"
        response.error_code = "CALIBRATION_ERROR"
        response.message = str(exc)

    return response.finalize()


def convert_ras_to_srh(payload: Dict[str, Any]) -> ToolResponse:
    """Convert RAS2D mesh/material inputs to SRH-2D format."""
    response = ToolResponse(status="running", message="Conversion started")

    try:
        req = validate_convert_request(payload)
        converter = pyHMT2D.Misc.RAS_to_SRH_Converter(
            req["ras_hdf_file"], req["terrain_tif_file"], req["srh_case_name"]
        )
        converter.convert_to_SRH()

        response.status = "succeeded"
        response.message = "Conversion completed"
        out_geom = Path(f"{req['srh_case_name']}.srhgeom")
        out_mat = Path(f"{req['srh_case_name']}.srhmat")
        response.artifacts.extend([_artifact(out_geom, with_hash=True), _artifact(out_mat, with_hash=True)])

        manifest = write_run_manifest(Path.cwd(), req, response)
        response.artifacts.append(_artifact(manifest, with_hash=True))
    except ValidationError as exc:
        response.status = "failed"
        response.error_code = "VALIDATION_ERROR"
        response.message = str(exc)
    except SystemExit as exc:
        response.status = "failed"
        response.error_code = "MODEL_EXIT"
        response.message = f"Model exited during conversion: {exc}"
    except Exception as exc:  # noqa: BLE001
        response.status = "failed"
        response.error_code = "CONVERSION_ERROR"
        response.message = str(exc)

    return response.finalize()


def summarize_artifacts(payload: Dict[str, Any]) -> ToolResponse:
    """Summarize matching artifacts under a directory tree."""
    response = ToolResponse(status="running", message="Artifact summary started")

    try:
        req = validate_result_summary_request(payload)
        root = Path(req["root"])
        suffixes = tuple(req["suffixes"])
        max_files = req["max_files"]

        matches: List[Path] = []
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            if path.suffix.lower() in suffixes:
                matches.append(path)
            if len(matches) >= max_files:
                break

        response.artifacts = [_artifact(p, with_hash=False) for p in matches]
        response.metrics = {
            "root": str(root),
            "suffixes": list(suffixes),
            "file_count": len(matches),
            "max_files": max_files,
        }
        response.status = "succeeded"
        response.message = "Artifact summary completed"

        manifest = write_run_manifest(root, req, response)
        response.artifacts.append(_artifact(manifest, with_hash=True))
    except ValidationError as exc:
        response.status = "failed"
        response.error_code = "VALIDATION_ERROR"
        response.message = str(exc)
    except Exception as exc:  # noqa: BLE001
        response.status = "failed"
        response.error_code = "SUMMARY_ERROR"
        response.message = str(exc)

    return response.finalize()


def run_single_simulation(payload: Dict[str, Any]) -> ToolResponse:
    """Validation-first run-model wrapper."""
    response = ToolResponse(status="running", message="Run-model request accepted")

    try:
        req = validate_run_model_request(payload)
        response.metrics["normalized_request"] = req
        response.status = "succeeded"
        response.message = (
            "Run-model request validated. Use existing model APIs to execute "
            "solver-specific runs in this environment."
        )

        manifest = write_run_manifest(Path(req["project_file"]).parent, req, response)
        response.artifacts.append(_artifact(manifest, with_hash=True))
    except ValidationError as exc:
        response.status = "failed"
        response.error_code = "VALIDATION_ERROR"
        response.message = str(exc)
    except JobExecutionError as exc:
        response.status = "failed"
        response.error_code = "JOB_EXECUTION_ERROR"
        response.message = str(exc)

    return response.finalize()


def modify_srh_case_parameters(payload: Dict[str, Any]) -> ToolResponse:
    """Modify Manning's n and inlet discharges in an existing SRH-2D case."""
    response = ToolResponse(status="running", message="SRH-2D case modification started")

    try:
        req = validate_modify_srh_case_request(payload)
        control_file = Path(req["srh_control_file"])

        data_obj = pyHMT2D.SRH_2D.SRH_2D_Data(str(control_file))

        if req["material_ids"]:
            data_obj.modify_ManningsNs(req["material_ids"], req["manning_n_values"], [])
        if req["inlet_bc_ids"]:
            data_obj.modify_InletQ(req["inlet_bc_ids"], req["inlet_q_values"])

        saved_control = _save_modified_srh_control(data_obj, req["output_control_file"])
        response.status = "succeeded"
        response.message = "SRH-2D case parameters updated"
        response.metrics = {
            "control_type": data_obj.control_type,
            "num_manning_updates": len(req["material_ids"]),
            "num_inlet_q_updates": len(req["inlet_bc_ids"]),
            "saved_control_file": str(saved_control),
        }
        response.artifacts.append(_artifact(saved_control, with_hash=True))

        manifest = write_run_manifest(saved_control.parent, req, response)
        response.artifacts.append(_artifact(manifest, with_hash=True))
    except ValidationError as exc:
        response.status = "failed"
        response.error_code = "VALIDATION_ERROR"
        response.message = str(exc)
    except Exception as exc:  # noqa: BLE001
        response.status = "failed"
        response.error_code = "MODIFY_SRH_CASE_ERROR"
        response.message = str(exc)

    return response.finalize()


def _build_srh_model(model_cfg: Dict[str, Any]):
    required = ["version", "srh_pre_path", "srh_path", "extra_dll_path"]
    missing = [k for k in required if k not in model_cfg]
    if missing:
        raise ValidationError(f"srh_model missing required keys for execution: {missing}")

    return pyHMT2D.SRH_2D.SRH_2D_Model(
        model_cfg["version"],
        model_cfg["srh_pre_path"],
        model_cfg["srh_path"],
        model_cfg["extra_dll_path"],
        faceless=bool(model_cfg.get("faceless", True)),
    )


def run_monte_carlo_srh(payload: Dict[str, Any]) -> ToolResponse:
    """Generate Monte Carlo samples for Manning's n and inlet discharge and optionally run SRH-2D."""
    response = ToolResponse(status="running", message="SRH Monte Carlo workflow started")

    try:
        req = validate_monte_carlo_srh_request(payload)
        base_control = Path(req["srh_control_file"]).resolve()
        base_case_dir = base_control.parent
        out_dir = Path(req["output_dir"]).resolve()
        out_dir.mkdir(parents=True, exist_ok=True)

        rng = random.Random(req["seed"])
        sample_records: List[Dict[str, Any]] = []

        for sample_idx in range(req["n_samples"]):
            sample_name = f"sample_{sample_idx:04d}"
            sample_dir = out_dir / sample_name
            if sample_dir.exists():
                shutil.rmtree(sample_dir)
            shutil.copytree(base_case_dir, sample_dir)

            sample_control = sample_dir / base_control.name
            sample_data = pyHMT2D.SRH_2D.SRH_2D_Data(str(sample_control))

            sampled_material_ids: List[int] = []
            sampled_manning: List[float] = []
            for m in req["manning_n_ranges"]:
                sampled_material_ids.append(m["material_id"])
                sampled_manning.append(rng.uniform(m["min"], m["max"]))

            sampled_inlet_ids: List[int] = []
            sampled_inlet_q: List[float] = []
            for b in req["inlet_q_ranges"]:
                sampled_inlet_ids.append(b["bc_id"])
                sampled_inlet_q.append(rng.uniform(b["min"], b["max"]))

            if sampled_material_ids:
                sample_data.modify_ManningsNs(sampled_material_ids, sampled_manning, [])
            if sampled_inlet_ids:
                sample_data.modify_InletQ(sampled_inlet_ids, sampled_inlet_q)

            _save_modified_srh_control(sample_data, None)

            run_status = "not-run"
            run_error = ""
            if req["run_simulations"]:
                with _pushd(sample_dir):
                    model = _build_srh_model(req["srh_model"])
                    try:
                        model.init_model()
                        model.open_project(sample_control.name)
                        model.run_pre_model()
                        model.run_model()
                        model.close_project()
                        model.exit_model()
                        run_status = "succeeded"
                    except Exception as exc:  # noqa: BLE001
                        run_status = "failed"
                        run_error = str(exc)

            sample_records.append(
                {
                    "sample_name": sample_name,
                    "sample_dir": str(sample_dir),
                    "material_ids": sampled_material_ids,
                    "manning_n_values": sampled_manning,
                    "inlet_bc_ids": sampled_inlet_ids,
                    "inlet_q_values": sampled_inlet_q,
                    "run_status": run_status,
                    "run_error": run_error,
                }
            )

        samples_json = out_dir / "monte_carlo_samples.json"
        with open(samples_json, "w", encoding="utf-8") as f:
            json.dump(sample_records, f, indent=2)

        samples_csv = out_dir / "monte_carlo_samples.csv"
        with open(samples_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
                    "sample_name",
                    "sample_dir",
                    "material_ids",
                    "manning_n_values",
                    "inlet_bc_ids",
                    "inlet_q_values",
                    "run_status",
                    "run_error",
                ]
            )
            for rec in sample_records:
                writer.writerow(
                    [
                        rec["sample_name"],
                        rec["sample_dir"],
                        rec["material_ids"],
                        rec["manning_n_values"],
                        rec["inlet_bc_ids"],
                        rec["inlet_q_values"],
                        rec["run_status"],
                        rec["run_error"],
                    ]
                )

        response.status = "succeeded"
        response.message = "SRH Monte Carlo workflow completed"
        response.metrics = {
            "n_samples": req["n_samples"],
            "output_dir": str(out_dir),
            "run_simulations": req["run_simulations"],
            "failed_runs": sum(1 for s in sample_records if s["run_status"] == "failed"),
        }
        response.artifacts.extend(
            [
                _artifact(samples_json, with_hash=True),
                _artifact(samples_csv, with_hash=True),
            ]
        )

        manifest = write_run_manifest(out_dir, req, response)
        response.artifacts.append(_artifact(manifest, with_hash=True))
    except ValidationError as exc:
        response.status = "failed"
        response.error_code = "VALIDATION_ERROR"
        response.message = str(exc)
    except Exception as exc:  # noqa: BLE001
        response.status = "failed"
        response.error_code = "MONTE_CARLO_SRH_ERROR"
        response.message = str(exc)

    return response.finalize()
