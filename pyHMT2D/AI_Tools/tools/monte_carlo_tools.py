# -*- coding: utf-8 -*-
"""
Group 7 — Monte Carlo simulation tools.

Based directly on the example scripts in
examples/Monte_Carlo/SRH_2D/ and examples/Monte_Carlo/RAS_2D/.

The MC workflow:
  1. generate_mc_samples  : draw N parameter samples from truncated-normal
                            (or uniform) distributions.
  2. run_monte_carlo      : copy the base case N times, apply one sample per
                            case, run all cases (parallel or serial), collect
                            last-timestep VTK results.
  3. get_mc_statistics    : compute exceedance probabilities at probe points
                            and (optionally) over the full spatial domain.

param_specs structure (extends the calibration schema with distribution info):
[
  {
    "type": "manning_n",           # "manning_n" | "inlet_q" | "exit_h"
    "material_id": 2,              # required for manning_n
    "material_name": "channel",    # optional
    "bc_id": None,                 # required for inlet_q / exit_h
    "distribution": "truncated_normal",   # or "uniform"
    "mean": 0.04,                  # used by truncated_normal
    "std":  0.005,                 # used by truncated_normal
    "min":  0.03,                  # lower bound (both distributions)
    "max":  0.05                   # upper bound (both distributions)
  },
  ...
]
"""

from __future__ import annotations

import contextlib
import csv
import json
import math
import multiprocessing
import os
import shutil
import sys as _sys
import time
from typing import List, Optional

import numpy as np

from pyHMT2D.AI_Tools.state import get_session


def _ok(data, message: str) -> dict:
    return {"status": "ok", "data": data, "message": message}


def _err(message: str) -> dict:
    return {"status": "error", "data": None, "message": message}


def _write_log(log_path: str, msg: str) -> None:
    """Write a timestamped line to stderr and append to the log file."""
    line = f"{time.strftime('%H:%M:%S')}  {msg}"
    print(line, file=_sys.stderr, flush=True)
    try:
        with open(log_path, "a") as fh:
            fh.write(line + "\n")
    except OSError:
        pass


class _TeeStream:
    """File-like object that writes to a log file and stderr simultaneously."""

    def __init__(self, log_path: str):
        self._log_path = log_path
        self._stderr = _sys.stderr

    def write(self, text: str) -> int:
        self._stderr.write(text)
        self._stderr.flush()
        try:
            with open(self._log_path, "a") as fh:
                fh.write(text)
        except OSError:
            pass
        return len(text)

    def flush(self) -> None:
        self._stderr.flush()

    def fileno(self):
        return self._stderr.fileno()


@contextlib.contextmanager
def _capture_model_output(log_path: str):
    """Redirect sys.stdout to log file + stderr for duration of block."""
    tee = _TeeStream(log_path)
    _saved = _sys.stdout
    _sys.stdout = tee
    try:
        yield
    finally:
        _sys.stdout = _saved


# ── helpers ───────────────────────────────────────────────────────────────────

# Fuzzy prefix patterns for auto-detecting WSE variables (case-insensitive).
_WSE_PREFIXES = ["water_elev", "water surface", "watersurface", "wse"]


def _get_vtk_variable_names(vtk_file: str) -> list:
    """Return all cell-data variable names in a VTK unstructured grid file."""
    from pyHMT2D.Misc import vtkHandler
    handler = vtkHandler()
    reader = handler.readVTK_UnstructuredGrid(vtk_file)
    ugrid = reader.GetOutput()
    cd = ugrid.GetCellData()
    return [cd.GetArrayName(i) for i in range(cd.GetNumberOfArrays())]


def _resolve_variable(vtk_file: str, variable) -> str:
    """Return a valid variable name from the VTK file.

    If *variable* is provided and exists, return it unchanged.
    If it is None, scan all VTK variables and return the first whose
    lowercase name starts with any entry in ``_WSE_PREFIXES`` (fuzzy match).
    Raises RuntimeError with available names if nothing matches.
    """
    available = _get_vtk_variable_names(vtk_file)
    if variable is not None:
        if variable in available:
            return variable
        raise RuntimeError(
            f"variable '{variable}' not found in VTK. "
            f"Available variables: {available}"
        )
    for name in available:
        if any(name.lower().startswith(p) for p in _WSE_PREFIXES):
            return name
    raise RuntimeError(
        f"Could not auto-detect a water surface elevation variable. "
        f"Available variables: {available}. Pass variable explicitly."
    )


def _sample_from_spec(spec: dict, n: int, rng: np.random.Generator) -> np.ndarray:
    """Draw n samples from the distribution described in spec."""
    dist = spec.get("distribution", "truncated_normal")
    lo, hi = float(spec["min"]), float(spec["max"])

    if dist == "uniform":
        return rng.uniform(lo, hi, size=n)

    # truncated_normal
    mean = float(spec.get("mean", (lo + hi) / 2))
    std  = float(spec.get("std",  (hi - lo) / 6))

    from scipy.stats import truncnorm
    a = (lo - mean) / std
    b = (hi - mean) / std
    return truncnorm.rvs(a, b, loc=mean, scale=std, size=n, random_state=int(rng.integers(1e9)))


def _run_single_case(args: tuple) -> dict:
    """Worker function: run one MC case. Called by multiprocessing.Pool."""
    (
        case_id, n_total, base_case_dir, output_dir, model_type, project_basename,
        param_specs, param_values,
        srh_path, srh_pre_path, extra_dll_path, srh_2d_version, hecras_version, faceless,
        delete_case_dir, log_path,
    ) = args

    # Build absolute case_dir here so os.chdir is always unambiguous
    case_dir = os.path.join(os.path.abspath(output_dir), f"case_{case_id:06d}")
    result = {"case_id": case_id, "success": False, "vtk_file": None, "error": None}

    param_str = ", ".join(
        f"{s.get('material_name', s.get('type', 'p'))}={v:.4f}"
        for s, v in zip(param_specs, param_values)
    )
    _write_log(log_path, f"Case {case_id}/{n_total}  started   params: {param_str}")
    t_case = time.time()

    # pyHMT2D model output log — written alongside the MC progress log
    model_log = os.path.join(os.path.abspath(output_dir), "pyHMT2D.log")

    try:
        import pyHMT2D

        # Create directory and copy base case files and subdirectories for this run
        os.makedirs(case_dir, exist_ok=True)
        for fn in os.listdir(base_case_dir):
            src = os.path.join(base_case_dir, fn)
            dst = os.path.join(case_dir, fn)
            if os.path.isfile(src):
                shutil.copy2(src, dst)
            elif os.path.isdir(src):
                shutil.copytree(src, dst)

        run_project_file = os.path.join(case_dir, project_basename)

        with _capture_model_output(model_log):
            if model_type == "SRH-2D":
                run_data = pyHMT2D.SRH_2D.SRH_2D_Data(run_project_file)
            else:
                # RAS_2D_Data expects the plan result HDF file (e.g. Muncie2D.p01.hdf),
                # not the .prj project file.  The session's hdf_filename gives the basename.
                from pyHMT2D.AI_Tools.state import get_session as _get_session
                _sess = _get_session()
                if _sess.data is not None and hasattr(_sess.data, 'hdf_filename'):
                    base_hdf = os.path.basename(_sess.data.hdf_filename)
                else:
                    import glob as _glob
                    plan_hdfs = sorted(_glob.glob(os.path.join(case_dir, "*.p[0-9]*.hdf")))
                    if not plan_hdfs:
                        raise RuntimeError(
                            "No plan result HDF file found in case directory. "
                            "Open the project with get_project_info first."
                        )
                    base_hdf = os.path.basename(plan_hdfs[0])
                plan_hdf_abs = os.path.abspath(os.path.join(case_dir, base_hdf))
                if not os.path.isfile(plan_hdf_abs):
                    raise RuntimeError(
                        f"Plan result HDF not found in case directory: {plan_hdf_abs}"
                    )
                run_data = pyHMT2D.RAS_2D.RAS_2D_Data(plan_hdf_abs)

            for spec, val in zip(param_specs, param_values):
                ptype = spec["type"]
                if ptype == "manning_n":
                    mat_id = int(spec["material_id"])
                    name   = spec.get("material_name", f"Material_{mat_id}")
                    if model_type == "SRH-2D":
                        run_data.modify_ManningsNs([mat_id], [float(val)], [name])
                    else:
                        run_data.modify_ManningsN([mat_id], [float(val)], [name])
                elif ptype == "inlet_q":
                    run_data.modify_InletQ([int(spec["bc_id"])], [float(val)])
                elif ptype == "exit_h":
                    run_data.modify_ExitH([int(spec["bc_id"])], [float(val)])

            if model_type == "SRH-2D":
                ctrl = (run_data.srhhydro_obj if run_data.control_type == "SRHHydro"
                        else run_data.srhsif_obj)
                ctrl.save_as()

            if model_type == "SRH-2D":
                model = pyHMT2D.SRH_2D.SRH_2D_Model(
                    srh_2d_version, srh_pre_path, srh_path, extra_dll_path or "", faceless=True
                )
                model.init_model()
                model.set_simulation_case(run_data)

                _orig_dir = os.getcwd()
                os.chdir(case_dir)
                try:
                    model.run_pre_model()
                    model.run_model(sleepTime=5.0, bShowProgress=False)
                finally:
                    os.chdir(_orig_dir)

                import glob as _glob
                xmdf_files = _glob.glob(os.path.join(case_dir, "*.h5"))
                if not xmdf_files:
                    raise RuntimeError("No XMDF result file found.")
                run_data.readSRHXMDFFile(xmdf_files[0], bNodal=False)
                run_data.outputXMDFDataToVTK(bNodal=False, lastTimeStep=True, dir=case_dir)

            else:  # HEC-RAS
                import glob as _glob
                run_prj_files = _glob.glob(os.path.join(case_dir, "*.prj"))
                run_prj_file = run_prj_files[0] if run_prj_files else run_project_file
                model = pyHMT2D.RAS_2D.HEC_RAS_Model(version=hecras_version, faceless=faceless)
                model.init_model()
                model.open_project(run_prj_file)
                model.run_model()
                rd = model.get_simulation_case()
                rd.saveHEC_RAS2D_results_to_VTK(
                    lastTimeStep=True, fileNameBase="ras_result", dir=case_dir
                )
                model.close_project()
                model.exit_model()

        import glob as _glob
        vtk_files = sorted(_glob.glob(os.path.join(case_dir, "*.vtk")))
        if not vtk_files:
            raise RuntimeError("Simulation produced no VTK output.")

        vtk_dest = case_dir.rstrip("/\\") + ".vtk"
        shutil.copy2(vtk_files[-1], vtk_dest)
        result["vtk_file"] = vtk_dest
        result["success"] = True
        elapsed = round(time.time() - t_case, 1)
        _write_log(log_path, f"Case {case_id}/{n_total}  completed  elapsed={elapsed}s")

    except Exception as exc:
        result["error"] = str(exc)
        elapsed = round(time.time() - t_case, 1)
        _write_log(log_path, f"Case {case_id}/{n_total}  FAILED ({elapsed}s): {exc}")

    finally:
        if delete_case_dir and os.path.isdir(case_dir):
            shutil.rmtree(case_dir, ignore_errors=True)

    return result


# ── tool 23 ───────────────────────────────────────────────────────────────────

def generate_mc_samples(
    param_specs: list,
    n_samples: int = 100,
    random_seed: int = 42,
    output_csv: Optional[str] = None,
) -> dict:
    """Generate Monte Carlo parameter samples from specified distributions.

    Uses truncated normal or uniform distributions. Saves samples to a CSV
    file for auditing and reproducing results.

    Parameters
    ----------
    param_specs : list of dict
        Each entry describes one parameter and its distribution (see module
        docstring for full schema). Required keys: "type", "min", "max".
        Optional: "distribution" (default "truncated_normal"), "mean", "std".
    n_samples : int
        Number of Monte Carlo samples to generate. Default 100.
    random_seed : int
        Random seed for reproducibility. Default 42.
    output_csv : str, optional
        Path to write the sample CSV. Defaults to 'mc_samples.csv' in the
        project directory (or current directory if no project is open).

    Returns
    -------
    dict with "data" containing
        {"n_samples": int, "param_names": [...], "samples": [[...], ...],
         "sample_csv": str}
        where samples[i] is a list of parameter values for the i-th run.
    """
    if n_samples < 1:
        return _err("n_samples must be at least 1.")

    rng = np.random.default_rng(random_seed)
    param_names = [
        f"{s['type']}_{s.get('material_id', s.get('bc_id', i))}"
        for i, s in enumerate(param_specs)
    ]

    try:
        cols = []
        for spec in param_specs:
            cols.append(_sample_from_spec(spec, n_samples, rng).tolist())

        # Transpose: samples[i] = values for run i
        samples = [[cols[j][i] for j in range(len(param_specs))] for i in range(n_samples)]

        # Determine output path
        session = get_session()
        if output_csv is None:
            base_dir = (
                os.path.dirname(os.path.abspath(session.project_file))
                if session.is_open() and session.project_file
                else os.getcwd()
            )
            output_csv = os.path.join(base_dir, "mc_samples.csv")

        with open(output_csv, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(param_names)
            for row in samples:
                writer.writerow([round(v, 6) for v in row])

        return _ok(
            {
                "n_samples": n_samples,
                "param_names": param_names,
                "samples": [[round(v, 6) for v in row] for row in samples],
                "sample_csv": output_csv,
            },
            f"Generated {n_samples} samples for {len(param_specs)} parameter(s). "
            f"Saved to {output_csv}.",
        )

    except Exception as exc:
        return _err(f"Sample generation failed: {exc}")


# ── tool 24 ───────────────────────────────────────────────────────────────────

def run_monte_carlo(
    base_case_dir: str,
    param_specs: list,
    n_samples: int = 100,
    n_processes: int = 1,
    random_seed: int = 42,
    sample_csv: Optional[str] = None,
    delete_cases: bool = True,
    output_dir: Optional[str] = None,
    srh_path: Optional[str] = None,
    srh_pre_path: Optional[str] = None,
    extra_dll_path: Optional[str] = None,
    srh_2d_version: str = "3.7",
    hecras_version: str = "6.6",
    faceless: bool = True,
) -> dict:
    """Run a Monte Carlo ensemble of hydraulic simulations.

    Creates N copies of the base case, applies one sample of parameters to each,
    runs all cases (in parallel with multiprocessing or serially), and collects
    last-timestep VTK results.

    Follows the workflow of examples/Monte_Carlo/SRH_2D/ and
    examples/Monte_Carlo/RAS_2D/.

    Parameters
    ----------
    base_case_dir : str
        Directory containing the template model files (the "base case").
        All files in this directory are copied for each MC run.
    param_specs : list of dict
        Parameter distributions (see module docstring for schema).
    n_samples : int
        Number of Monte Carlo samples (model runs). Default 100.
    n_processes : int
        Number of parallel worker processes. Use 1 for serial execution.
        Default 1 (safe; set to number of CPU cores for speed).
    random_seed : int
        Random seed for reproducibility. Used only when ``sample_csv`` is not
        provided. Default 42.
    sample_csv : str, optional
        Path to a sample CSV previously written by ``generate_mc_samples()``.
        If provided, samples are read from this file instead of being drawn
        internally, and ``n_samples`` / ``random_seed`` are ignored.
        This ensures the runs use exactly the same samples that were previewed
        and audited in the generate step.
    delete_cases : bool
        If True, delete each case directory after copying its VTK result.
        Saves disk space. Default True.
    output_dir : str, optional
        Directory to store case directories and results. Defaults to a
        'mc_runs' subdirectory inside base_case_dir.
    srh_path, srh_pre_path, extra_dll_path : str
        SRH-2D solver paths. Required for SRH-2D runs.
    srh_2d_version : str
        SRH-2D version. Default "3.7".
    hecras_version : str
        HEC-RAS version. Default "6.6".

    Returns
    -------
    dict with "data" containing
        {"results_json": str, "n_successful": int, "n_failed": int,
         "vtk_files": [...], "case_dir": str, "elapsed_seconds": float}
    """
    session = get_session()
    try:
        session.require_open()
    except RuntimeError as e:
        return _err(str(e))

    # Use session's configured HEC-RAS version if not explicitly provided
    hecras_version = getattr(session, "_hecras_version", None) or hecras_version

    if not os.path.isdir(base_case_dir):
        return _err(f"base_case_dir not found: {base_case_dir}")

    if output_dir is None:
        output_dir = os.path.join(os.path.dirname(os.path.abspath(base_case_dir)), "mc_runs")
    os.makedirs(output_dir, exist_ok=True)

    # Progress log — overwrite any existing log from a previous run
    log_path = os.path.join(os.path.abspath(output_dir), "mc_progress.log")
    try:
        open(log_path, "w").close()
    except OSError:
        pass

    # Load or generate samples
    if sample_csv is not None:
        if not os.path.isfile(sample_csv):
            return _err(f"sample_csv not found: {sample_csv}")
        try:
            with open(sample_csv, "r") as fh:
                reader = csv.reader(fh)
                next(reader)  # skip header row
                samples = [[float(v) for v in row] for row in reader if row]
        except Exception as exc:
            return _err(f"Failed to read sample_csv: {exc}")
        n_samples = len(samples)
        _write_log(log_path, f"Loaded {n_samples} samples from {sample_csv}")
    else:
        rng = np.random.default_rng(random_seed)
        param_cols = [_sample_from_spec(spec, n_samples, rng) for spec in param_specs]
        samples = [[float(param_cols[j][i]) for j in range(len(param_specs))] for i in range(n_samples)]
        _write_log(log_path, f"Generated {n_samples} samples internally (seed={random_seed})")

    _write_log(log_path, f"MC started: {n_samples} runs, {n_processes} process(es)")

    project_basename = os.path.basename(session.project_file)
    abs_base_case_dir = os.path.abspath(base_case_dir)
    abs_output_dir = os.path.abspath(output_dir)
    t0 = time.time()

    # Build arguments list for workers — case directories are created inside
    # _run_single_case so they appear one at a time as each run starts.
    worker_args = []
    for i, sample in enumerate(samples):
        worker_args.append((
            i + 1, n_samples, abs_base_case_dir, abs_output_dir, session.model_type, project_basename,
            param_specs, sample,
            srh_path or getattr(session, "_srh_path", ""),
            srh_pre_path or getattr(session, "_srh_pre_path", ""),
            extra_dll_path or getattr(session, "_extra_dll_path", ""),
            srh_2d_version, hecras_version, faceless,
            delete_cases, log_path,
        ))

    # Run
    all_results = []
    if n_processes > 1:
        try:
            with multiprocessing.Pool(processes=n_processes, maxtasksperchild=1) as pool:
                all_results = pool.map(_run_single_case, worker_args)
        except Exception as exc:
            return _err(f"Parallel execution failed: {exc}. Try n_processes=1.")
    else:
        for args in worker_args:
            all_results.append(_run_single_case(args))

    # Summarise
    successful = [r for r in all_results if r["success"]]
    failed     = [r for r in all_results if not r["success"]]
    vtk_files  = [r["vtk_file"] for r in successful if r["vtk_file"]]

    elapsed = round(time.time() - t0, 1)
    _write_log(log_path, f"MC complete: {len(successful)}/{n_samples} succeeded in {elapsed}s")

    summary = {
        "n_samples": n_samples,
        "n_successful": len(successful),
        "n_failed": len(failed),
        "all_run_results": [r["case_id"] if r["success"] else -r["case_id"] for r in all_results],
        "successful_runs": [r["case_id"] for r in successful],
        "failed_runs":     [r["case_id"] for r in failed],
        "failed_errors":   {r["case_id"]: r["error"] for r in failed},
        "vtk_files":       vtk_files,
        "elapsed_seconds": elapsed,
    }

    # Save JSON
    results_json = os.path.join(output_dir, f"mc_results_{time.strftime('%Y%m%d_%H%M%S')}.json")
    with open(results_json, "w") as fh:
        json.dump(summary, fh, indent=2)

    return _ok(
        {
            "results_json": results_json,
            "n_successful": len(successful),
            "n_failed": len(failed),
            "vtk_files": vtk_files,
            "case_dir": output_dir,
            "elapsed_seconds": elapsed,
            "progress_log": log_path,
        },
        f"MC complete: {len(successful)}/{n_samples} runs succeeded in {elapsed}s. "
        f"Results in {output_dir}.",
    )


# ── tool 25 ───────────────────────────────────────────────────────────────────

def get_mc_statistics(
    results_json: str,
    variable: Optional[str] = None,
    observation_points: Optional[List[dict]] = None,
    exceedance_probabilities: Optional[List[float]] = None,
    output_vtk: Optional[str] = None,
) -> dict:
    """Compute exceedance-probability statistics from a Monte Carlo ensemble.

    Aggregates last-timestep VTK results across all successful runs, computes
    user-specified exceedance probabilities at probe points and (optionally)
    over the full spatial domain.

    Follows the post-processing pattern in
    examples/Monte_Carlo/SRH_2D/demo_SRH_2D_Monte_Carlo.py.

    Parameters
    ----------
    results_json : str
        Path to the JSON summary file written by run_monte_carlo().
    variable : str, optional
        VTK variable to analyse. If omitted, auto-detected from the VTK file
        using fuzzy prefix matching (water_elev*, wse*, water surface*).
    observation_points : list of dict, optional
        Point locations for CDF extraction. Each entry:
          {"name": str, "x": float, "y": float}
        If omitted, only spatial statistics are computed.
    exceedance_probabilities : list of float, optional
        Exceedance probability percentiles (0–100) to compute.
        Default [99, 90, 50, 10, 1].
    output_vtk : str, optional
        Path for the spatial-statistics VTK output file. Defaults to
        'mc_exceedance.vtk' in the same directory as results_json.

    Returns
    -------
    dict with "data" containing
        {"point_stats": [...], "spatial_vtk": str, "n_runs": int,
         "exceedance_probabilities": [...]}
    where point_stats[i] = {"name": str, "x": float, "y": float,
                             "values_at_probs": {prob: value, ...}}
    """
    if not os.path.isfile(results_json):
        return _err(f"Results JSON not found: {results_json}")

    if exceedance_probabilities is None:
        exceedance_probabilities = [99, 90, 50, 10, 1]

    _model_log = os.path.join(os.path.dirname(os.path.abspath(results_json)), "pyHMT2D.log")

    try:
        with open(results_json) as fh:
            summary = json.load(fh)

        vtk_files = summary.get("vtk_files", [])
        if not vtk_files:
            return _err("No VTK files listed in results JSON. Were any runs successful?")

        vtk_files = [f for f in vtk_files if os.path.isfile(f)]
        if not vtk_files:
            return _err("VTK files listed in results JSON no longer exist on disk.")

        point_stats = []
        spatial_vtk_path = output_vtk

        with _capture_model_output(_model_log):
            from pyHMT2D.Misc import vtkHandler
            import vtk as _vtk

            handler = vtkHandler()

            # Auto-detect variable name from first VTK file if not provided
            variable = _resolve_variable(vtk_files[0], variable)

            # --- Point-based statistics ---
            if observation_points:
                # Probe each VTK file at each point
                n_pts = len(observation_points)
                n_runs = len(vtk_files)
                probe_matrix = np.full((n_pts, n_runs), np.nan)

                for run_i, vf in enumerate(vtk_files):
                    try:
                        reader = handler.readVTK_UnstructuredGrid(vf)
                        xs = [float(p["x"]) for p in observation_points]
                        ys = [float(p["y"]) for p in observation_points]
                        pts_vtk = _vtk.vtkPoints()
                        for x, y in zip(xs, ys):
                            pts_vtk.InsertNextPoint(x, y, 0.0)
                        probed = handler.probeUnstructuredGridVTKOnPoints(
                            pts_vtk, reader, variable
                        )
                        if probed is not None and len(probed) >= 2:
                            for pi, v in enumerate(probed[1]):
                                probe_matrix[pi, run_i] = float(v)
                    except Exception:
                        pass

                for pi, pt in enumerate(observation_points):
                    vals = np.sort(probe_matrix[pi, np.isfinite(probe_matrix[pi])])[::-1]
                    n_valid = len(vals)
                    probs_at = {}
                    for prob in exceedance_probabilities:
                        idx = min(int(prob * (n_valid + 1) / 100), n_valid - 1) if n_valid > 0 else 0
                        probs_at[f"p{prob}"] = round(float(vals[idx]), 4) if n_valid > 0 else None
                    point_stats.append({
                        "name": pt.get("name", f"pt_{pi}"),
                        "x": pt["x"],
                        "y": pt["y"],
                        "n_valid_runs": n_valid,
                        "values_at_exceedance_probs": probs_at,
                    })

            # --- Spatial statistics (cell-by-cell) ---
            if spatial_vtk_path is None:
                spatial_vtk_path = os.path.join(
                    os.path.dirname(results_json), "mc_exceedance.vtk"
                )

            # Read first VTK to get mesh structure
            reader0 = handler.readVTK_UnstructuredGrid(vtk_files[0])
            ugrid0 = reader0.GetOutput()
            n_cells = ugrid0.GetNumberOfCells()

            # Collect variable values across all runs
            all_vals = np.full((n_cells, len(vtk_files)), np.nan)
            for run_i, vf in enumerate(vtk_files):
                try:
                    reader = handler.readVTK_UnstructuredGrid(vf)
                    arr = handler.get_uGRid_cell_field_with_name(reader, variable)
                    if arr is not None:
                        all_vals[:, run_i] = np.array(arr).flatten()[:n_cells]
                except Exception:
                    pass

            # Compute exceedance probabilities at each cell
            from vtk.util import numpy_support as VN
            import vtk as _vtk2

            out_ugrid = _vtk2.vtkUnstructuredGrid()
            out_ugrid.DeepCopy(ugrid0)

            for prob in exceedance_probabilities:
                cell_probs = np.full(n_cells, np.nan)
                for ci in range(n_cells):
                    run_vals = np.sort(all_vals[ci, np.isfinite(all_vals[ci])])[::-1]
                    n_valid = len(run_vals)
                    if n_valid > 0:
                        idx = min(int(prob * (n_valid + 1) / 100), n_valid - 1)
                        cell_probs[ci] = run_vals[idx]

                vtk_arr = VN.numpy_to_vtk(
                    np.nan_to_num(cell_probs, nan=0.0), deep=True
                )
                vtk_arr.SetName(f"{variable}_p{prob}")
                out_ugrid.GetCellData().AddArray(vtk_arr)

            writer = _vtk2.vtkUnstructuredGridWriter()
            writer.SetFileName(spatial_vtk_path)
            writer.SetInputData(out_ugrid)
            writer.Write()

        output_vtk = spatial_vtk_path

        return _ok(
            {
                "n_runs": len(vtk_files),
                "variable": variable,
                "exceedance_probabilities": exceedance_probabilities,
                "point_stats": point_stats,
                "spatial_vtk": output_vtk,
            },
            f"MC statistics computed for {len(vtk_files)} runs. "
            f"Spatial exceedance VTK: {output_vtk}.",
        )

    except Exception as exc:
        return _err(f"MC statistics failed: {exc}")
