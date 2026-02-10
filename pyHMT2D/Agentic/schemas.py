"""Validation utilities for agent-facing requests.

The validators intentionally avoid heavy optional dependencies to keep
pyHMT2D install footprint stable.
"""

from pathlib import Path
from typing import Any, Dict, List

from .errors import ValidationError

SUPPORTED_MODELS = {"Backwater-1D", "SRH-2D", "HEC-RAS"}


def _assert_keys(payload: Dict[str, Any], required_keys: List[str], name: str) -> None:
    missing = [k for k in required_keys if k not in payload]
    if missing:
        raise ValidationError(f"{name} missing required keys: {missing}")


def _validate_id_value_pairs(
    payload: Dict[str, Any], ids_key: str, values_key: str, name: str
) -> Dict[str, List[Any]]:
    ids = payload.get(ids_key, [])
    values = payload.get(values_key, [])

    if len(ids) != len(values):
        raise ValidationError(f"{name}: '{ids_key}' and '{values_key}' must have same length")

    normalized_ids: List[int] = []
    normalized_values: List[float] = []

    for idx, (id_v, value_v) in enumerate(zip(ids, values)):
        try:
            normalized_ids.append(int(id_v))
        except Exception as exc:  # noqa: BLE001
            raise ValidationError(f"{name}: invalid integer at {ids_key}[{idx}]={id_v}") from exc

        try:
            normalized_values.append(float(value_v))
        except Exception as exc:  # noqa: BLE001
            raise ValidationError(f"{name}: invalid float at {values_key}[{idx}]={value_v}") from exc

    return {ids_key: normalized_ids, values_key: normalized_values}


def validate_run_model_request(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Validate run-model request payload."""
    _assert_keys(payload, ["model", "project_file"], "RunModelRequest")

    model = payload["model"]
    if model not in SUPPORTED_MODELS:
        raise ValidationError(f"Unsupported model '{model}'. Supported: {sorted(SUPPORTED_MODELS)}")

    project_file = payload["project_file"]
    if not isinstance(project_file, str) or not project_file.strip():
        raise ValidationError("project_file must be a non-empty string")

    path = Path(project_file)
    if payload.get("require_existing_project", True) and not path.exists():
        raise ValidationError(f"project_file does not exist: {project_file}")

    normalized = {
        "model": model,
        "project_file": str(path),
        "require_existing_project": bool(payload.get("require_existing_project", True)),
        "metadata": payload.get("metadata", {}),
    }
    return normalized


def validate_calibration_request(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Validate calibration request payload."""
    _assert_keys(payload, ["config_json"], "CalibrationRequest")

    config_json = payload["config_json"]
    if not isinstance(config_json, str) or not config_json.strip():
        raise ValidationError("config_json must be a non-empty string")

    path = Path(config_json)
    if payload.get("require_existing_config", True) and not path.exists():
        raise ValidationError(f"config_json does not exist: {config_json}")

    return {
        "config_json": str(path),
        "require_existing_config": bool(payload.get("require_existing_config", True)),
        "dry_run": bool(payload.get("dry_run", False)),
        "metadata": payload.get("metadata", {}),
    }


def validate_convert_request(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Validate RAS-to-SRH conversion request payload."""
    _assert_keys(payload, ["ras_hdf_file", "terrain_tif_file", "srh_case_name"], "ConvertRequest")

    ras_file = Path(payload["ras_hdf_file"])
    terrain_file = Path(payload["terrain_tif_file"])
    case_name = payload["srh_case_name"]

    if not case_name or not isinstance(case_name, str):
        raise ValidationError("srh_case_name must be a non-empty string")
    if payload.get("require_existing_inputs", True):
        if not ras_file.exists():
            raise ValidationError(f"ras_hdf_file does not exist: {ras_file}")
        if not terrain_file.exists():
            raise ValidationError(f"terrain_tif_file does not exist: {terrain_file}")

    return {
        "ras_hdf_file": str(ras_file),
        "terrain_tif_file": str(terrain_file),
        "srh_case_name": case_name,
        "require_existing_inputs": bool(payload.get("require_existing_inputs", True)),
        "metadata": payload.get("metadata", {}),
    }


def validate_result_summary_request(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Validate artifact/result summary request payload."""
    _assert_keys(payload, ["root"], "ResultSummaryRequest")

    root = Path(payload["root"])
    if payload.get("require_existing_root", True) and not root.exists():
        raise ValidationError(f"root path does not exist: {root}")

    suffixes = payload.get("suffixes", [".json", ".hdf", ".vtk", ".srhgeom", ".srhmat", ".log"])
    if not isinstance(suffixes, list) or not all(isinstance(s, str) for s in suffixes):
        raise ValidationError("suffixes must be a list of strings")

    return {
        "root": str(root),
        "suffixes": suffixes,
        "require_existing_root": bool(payload.get("require_existing_root", True)),
        "max_files": int(payload.get("max_files", 500)),
    }


def validate_modify_srh_case_request(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Validate request for modifying SRH-2D case parameters."""
    _assert_keys(payload, ["srh_control_file"], "ModifySRHCaseRequest")

    control_file = Path(payload["srh_control_file"])
    if payload.get("require_existing_control", True) and not control_file.exists():
        raise ValidationError(f"srh_control_file does not exist: {control_file}")

    material_pairs = _validate_id_value_pairs(
        payload, "material_ids", "manning_n_values", "ModifySRHCaseRequest"
    )
    inlet_pairs = _validate_id_value_pairs(
        payload, "inlet_bc_ids", "inlet_q_values", "ModifySRHCaseRequest"
    )

    output_control_file = payload.get("output_control_file")
    if output_control_file is not None and not isinstance(output_control_file, str):
        raise ValidationError("output_control_file must be a string if provided")

    return {
        "srh_control_file": str(control_file),
        "material_ids": material_pairs["material_ids"],
        "manning_n_values": material_pairs["manning_n_values"],
        "inlet_bc_ids": inlet_pairs["inlet_bc_ids"],
        "inlet_q_values": inlet_pairs["inlet_q_values"],
        "output_control_file": output_control_file,
        "metadata": payload.get("metadata", {}),
    }


def validate_monte_carlo_srh_request(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Validate request for Monte Carlo sample generation and SRH-2D runs."""
    _assert_keys(payload, ["srh_control_file", "n_samples", "output_dir"], "MonteCarloSRHRequest")

    control_file = Path(payload["srh_control_file"])
    if payload.get("require_existing_control", True) and not control_file.exists():
        raise ValidationError(f"srh_control_file does not exist: {control_file}")

    try:
        n_samples = int(payload["n_samples"])
    except Exception as exc:  # noqa: BLE001
        raise ValidationError("n_samples must be an integer") from exc
    if n_samples <= 0:
        raise ValidationError("n_samples must be > 0")

    output_dir = Path(payload["output_dir"])

    manning_ranges = payload.get("manning_n_ranges", [])
    inlet_ranges = payload.get("inlet_q_ranges", [])

    def _normalize_ranges(items: List[Dict[str, Any]], id_key: str, name: str) -> List[Dict[str, Any]]:
        normalized = []
        for i, item in enumerate(items):
            if not isinstance(item, dict):
                raise ValidationError(f"{name}[{i}] must be an object")
            if id_key not in item or "min" not in item or "max" not in item:
                raise ValidationError(f"{name}[{i}] must contain {id_key}, min, max")

            try:
                item_id = int(item[id_key])
                min_v = float(item["min"])
                max_v = float(item["max"])
            except Exception as exc:  # noqa: BLE001
                raise ValidationError(f"{name}[{i}] has invalid numeric value") from exc

            if min_v > max_v:
                raise ValidationError(f"{name}[{i}] min cannot be larger than max")

            normalized.append({id_key: item_id, "min": min_v, "max": max_v})
        return normalized

    normalized_manning = _normalize_ranges(manning_ranges, "material_id", "manning_n_ranges")
    normalized_inlet = _normalize_ranges(inlet_ranges, "bc_id", "inlet_q_ranges")

    return {
        "srh_control_file": str(control_file),
        "n_samples": n_samples,
        "output_dir": str(output_dir),
        "manning_n_ranges": normalized_manning,
        "inlet_q_ranges": normalized_inlet,
        "seed": int(payload.get("seed", 123)),
        "run_simulations": bool(payload.get("run_simulations", False)),
        "srh_model": payload.get("srh_model", {}),
        "metadata": payload.get("metadata", {}),
    }
