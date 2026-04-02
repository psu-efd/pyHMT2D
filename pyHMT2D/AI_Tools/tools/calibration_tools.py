# -*- coding: utf-8 -*-
"""
Group 6 — Calibration tools.

Based directly on the example scripts in
examples/calibration/SRH-2D/ and examples/calibration/RAS-2D/.

The calibration loop pattern:
  1. Load observations from a CSV file.
  2. For each candidate parameter set: copy base case → modify params →
     run model → export VTK → probe VTK at observation points → RMSE.
  3. Optimizer (GP-minimize or Nelder-Mead) drives the loop.

Three tools are exposed:
  - load_observations   : validate and preview the observation CSV.
  - evaluate_parameters : single forward run with given parameter values.
  - run_calibration     : full automated optimization loop.

param_specs structure (same schema used by Monte Carlo tools):
[
  {
    "type": "manning_n",          # required: "manning_n" | "inlet_q" | "exit_h"
    "material_id": 2,             # required for manning_n
    "material_name": "channel",   # optional, for logging
    "bc_id": None,                # required for inlet_q / exit_h
    "initial": 0.04,              # starting value
    "min": 0.02,                  # lower bound
    "max": 0.08                   # upper bound
  },
  ...
]

observation CSV format (matching examples/calibration/*/HWMs.dat):
  # comment lines start with #
  x, y, observed_value          (no header; 3 columns = scalar WSE/depth)
  OR
  name, x, y, observed_value    (4 columns)
"""

from __future__ import annotations

import contextlib
import copy
import csv
import math
import os
import shutil
import sys as _sys
import time
from typing import List, Optional

import numpy as np

from pyHMT2D.AI_Tools.state import get_session


class _TeeStream:
    """File-like object that writes to a log file and to stderr simultaneously."""

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
    """Context manager: redirect sys.stdout to log file + stderr for duration."""
    tee = _TeeStream(log_path)
    _saved = _sys.stdout
    _sys.stdout = tee
    try:
        yield
    finally:
        _sys.stdout = _saved


def _ok(data, message: str) -> dict:
    return {"status": "ok", "data": data, "message": message}


def _err(message: str) -> dict:
    return {"status": "error", "data": None, "message": message}


# ── helpers ───────────────────────────────────────────────────────────────────

def _parse_observation_csv(csv_file: str) -> tuple:
    """Return (names, xs, ys, values) lists from the observation CSV."""
    names, xs, ys, vals = [], [], [], []
    with open(csv_file, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) == 3:
                names.append(f"pt_{len(names)+1}")
                xs.append(float(parts[0]))
                ys.append(float(parts[1]))
                vals.append(float(parts[2]))
            elif len(parts) >= 4:
                names.append(parts[0])
                xs.append(float(parts[1]))
                ys.append(float(parts[2]))
                vals.append(float(parts[3]))
    return names, xs, ys, vals


def _apply_params(session_data, model_type: str, param_specs: list, values: list) -> None:
    """Apply a vector of parameter values to the open model data object."""
    for spec, val in zip(param_specs, values):
        ptype = spec["type"]
        if ptype == "manning_n":
            mat_id = int(spec["material_id"])
            mat_name = spec.get("material_name", f"Material_{mat_id}")
            if model_type == "SRH-2D":
                session_data.modify_ManningsNs([mat_id], [float(val)], [mat_name])
            else:
                session_data.modify_ManningsN([mat_id], [float(val)], [mat_name])
        elif ptype == "inlet_q":
            bc_id = int(spec["bc_id"])
            session_data.modify_InletQ([bc_id], [float(val)])
        elif ptype == "exit_h":
            bc_id = int(spec["bc_id"])
            session_data.modify_ExitH([bc_id], [float(val)])


# Fuzzy prefix patterns for auto-detecting WSE variables (case-insensitive).
# A VTK variable matches if its lowercase form starts with any of these.
_WSE_PREFIXES = ["water_elev", "water surface", "watersurface", "wse"]


def _get_vtk_variable_names(vtk_file: str) -> list:
    """Return all cell-data variable names in a VTK unstructured grid file."""
    from pyHMT2D.Misc import vtkHandler
    handler = vtkHandler()
    reader = handler.readVTK_UnstructuredGrid(vtk_file)
    ugrid = reader.GetOutput()
    cd = ugrid.GetCellData()
    return [cd.GetArrayName(i) for i in range(cd.GetNumberOfArrays())]


def _resolve_result_variable(vtk_file: str, result_variable) -> str:
    """Return a valid variable name from the VTK file.

    If *result_variable* is provided and exists, return it unchanged.
    If it is None, scan all VTK variables and return the first whose
    lowercase name starts with any entry in ``_WSE_PREFIXES`` (fuzzy match).
    Raises RuntimeError with available names if nothing matches.
    """
    available = _get_vtk_variable_names(vtk_file)
    if result_variable is not None:
        if result_variable in available:
            return result_variable
        raise RuntimeError(
            f"result_variable '{result_variable}' not found in VTK. "
            f"Available variables: {available}"
        )
    for name in available:
        low = name.lower()
        if any(low.startswith(p) for p in _WSE_PREFIXES):
            return name
    raise RuntimeError(
        f"Could not auto-detect a water surface elevation variable. "
        f"Available variables: {available}. "
        f"Pass result_variable explicitly."
    )


def _probe_wse_at_points(vtk_file: str, xs: list, ys: list, var_name: str) -> list:
    """Probe a VTK file at (x, y) points. Returns list of floated values."""
    import vtk as _vtk
    from pyHMT2D.Misc import vtkHandler

    handler = vtkHandler()
    reader = handler.readVTK_UnstructuredGrid(vtk_file)

    pts = _vtk.vtkPoints()
    for xi, yi in zip(xs, ys):
        pts.InsertNextPoint(float(xi), float(yi), 0.0)

    probed = handler.probeUnstructuredGridVTKOnPoints(pts, reader, var_name)
    if probed is None or len(probed) < 2:
        raise RuntimeError(f"Probing failed for variable '{var_name}'.")
    return [float(v) for v in probed[1]]


def _rmse(simulated: list, measured: list) -> float:
    diffs = [(s - m) ** 2 for s, m in zip(simulated, measured) if math.isfinite(s)]
    return math.sqrt(sum(diffs) / len(diffs)) if diffs else float("nan")


# ── tool 20 ───────────────────────────────────────────────────────────────────

def check_observation_format(csv_file: str) -> dict:
    """Validate the format of a high-water-mark (HWM) observation file.

    Expected format
    ---------------
    - Lines beginning with ``#`` are treated as comments and ignored.
    - Every data row must have exactly **3 columns**: ``x, y, wse``
      OR exactly **4 columns**: ``name, x, y, wse``.
    - Mixed row widths (some rows 3-column, others 4-column) are not allowed.
    - ``x`` and ``y`` must be finite floats (coordinates in the model CRS).
    - ``wse`` (water surface elevation) must be a finite float.
    - The file must contain at least one data row.

    Parameters
    ----------
    csv_file : str
        Path to the observation file (e.g. ``HWMs.dat``).

    Returns
    -------
    dict with ``"data"`` containing
        ``{"n_points": int, "columns": str, "preview": [...]}``
        where ``columns`` is ``"x,y,wse"`` or ``"name,x,y,wse"`` and
        ``preview`` shows the first (up to 5) parsed rows.
    """
    if not os.path.isfile(csv_file):
        return _err(f"File not found: {csv_file}")

    errors = []
    rows = []
    expected_ncols = None

    try:
        with open(csv_file, "r") as fh:
            for lineno, raw in enumerate(fh, start=1):
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue

                parts = [p.strip() for p in line.split(",")]
                ncols = len(parts)

                # Determine / enforce consistent column count
                if expected_ncols is None:
                    if ncols not in (3, 4):
                        errors.append(
                            f"Line {lineno}: expected 3 or 4 columns, got {ncols}: {line!r}"
                        )
                        continue
                    expected_ncols = ncols
                elif ncols != expected_ncols:
                    errors.append(
                        f"Line {lineno}: inconsistent column count "
                        f"(expected {expected_ncols}, got {ncols}): {line!r}"
                    )
                    continue

                # Parse numeric fields
                try:
                    if ncols == 3:
                        x, y, wse = float(parts[0]), float(parts[1]), float(parts[2])
                        name = f"pt_{len(rows)+1}"
                    else:
                        name = parts[0]
                        x, y, wse = float(parts[1]), float(parts[2]), float(parts[3])
                except ValueError as ve:
                    errors.append(f"Line {lineno}: non-numeric value — {ve}: {line!r}")
                    continue

                # Check finite values
                if not (math.isfinite(x) and math.isfinite(y)):
                    errors.append(f"Line {lineno}: x or y is not finite ({x}, {y})")
                    continue
                if not math.isfinite(wse):
                    errors.append(f"Line {lineno}: wse is not finite ({wse})")
                    continue

                rows.append({"name": name, "x": x, "y": y, "wse": wse})

    except Exception as exc:
        return _err(f"Failed to read file: {exc}")

    if not rows and not errors:
        return _err("File contains no data rows (only comments or blank lines).")

    if errors:
        return _err(
            f"Format errors found in {csv_file} "
            f"({len(errors)} issue(s)):\n" + "\n".join(errors)
        )

    col_label = "x,y,wse" if expected_ncols == 3 else "name,x,y,wse"
    return _ok(
        {
            "n_points": len(rows),
            "columns": col_label,
            "preview": rows[:5],
        },
        f"Format OK — {len(rows)} HWM point(s), columns: {col_label}.",
    )


# ── tool 21 ───────────────────────────────────────────────────────────────────

def load_observations(csv_file: str) -> dict:
    """Load and validate an observation CSV file (high water marks, WSE, etc.).

    Expected format (matching pyHMT2D calibration examples):
      - Lines starting with # are comments and are ignored.
      - 3-column form:  x, y, observed_value
      - 4-column form:  name, x, y, observed_value

    Parameters
    ----------
    csv_file : str
        Path to the CSV file containing field observations.

    Returns
    -------
    dict with "data" containing a list of
        {"name": str, "x": float, "y": float, "observed_value": float}
    """
    if not os.path.isfile(csv_file):
        return _err(f"Observation file not found: {csv_file}")

    try:
        names, xs, ys, vals = _parse_observation_csv(csv_file)
        if not names:
            return _err(f"No valid data rows found in {csv_file}.")

        records = [
            {"name": n, "x": x, "y": y, "observed_value": v}
            for n, x, y, v in zip(names, xs, ys, vals)
        ]
        return _ok(
            records,
            f"Loaded {len(records)} observation point(s) from {csv_file}.",
        )
    except Exception as exc:
        return _err(f"Failed to parse observation file: {exc}")


# ── tool 21 ───────────────────────────────────────────────────────────────────

def evaluate_parameters(
    param_specs: list,
    observation_csv: str,
    result_variable: Optional[str] = None,
    base_case_dir: Optional[str] = None,
    work_dir: Optional[str] = None,
    run_dir: Optional[str] = None,
    srh_path: Optional[str] = None,
    srh_pre_path: Optional[str] = None,
    extra_dll_path: Optional[str] = None,
    version: str = "3.7",
    hecras_version: str = "6.6",
) -> dict:
    """Run a single forward simulation with the given parameters and return RMSE.

    This is the inner loop of any calibration strategy. It:
      1. Copies the base case files into a local run subdirectory.
      2. Applies the parameter values to the copied model data.
      3. Saves the modified input files.
      4. Executes the simulation.
      5. Exports results to VTK.
      6. Probes the VTK at observation point locations.
      7. Computes and returns RMSE against observed values.

    Use run_calibration() for automated optimization, or call this tool
    directly when you want to control the optimization loop yourself.

    Parameters
    ----------
    param_specs : list of dict
        Each entry: {"type": "manning_n", "material_id": int,
                     "material_name": str, "value": float}.
        The "value" key holds the value to test in this call.
    observation_csv : str
        Path to CSV with observed values (see load_observations()).
    result_variable : str
        VTK variable name to compare against observations. Default "Water_Elev_m".
    base_case_dir : str, optional
        Directory containing the original (unmodified) model files to copy
        from for each run. Defaults to the directory of the open project file.
        Use this to keep a clean ``base_case/`` subdirectory separate from
        your working directory.
    work_dir : str, optional
        Parent directory in which run subdirectories (``calib_run_XXXXXXXX``)
        are created. Defaults to the current working directory. Using a local
        path (rather than the system temp folder) avoids path issues with
        SRH-2D Pre when the project uses relative file references.
    run_dir : str, optional
        Explicit directory to use for this run, overriding ``work_dir``.
        If omitted a unique subdirectory is created inside ``work_dir``.
    srh_path, srh_pre_path, extra_dll_path : str, optional
        SRH-2D executable paths (inherited from session if not provided).
    version : str
        SRH-2D version. Default "3.3".
    hecras_version : str
        HEC-RAS version. Default "6.6".

    Returns
    -------
    dict with "data" containing
        {"rmse": float, "per_point": [{"name", "x", "y", "simulated", "measured", "error"}],
         "param_values": {name: value}, "elapsed_seconds": float, "run_dir": str}
    """
    session = get_session()
    try:
        session.require_open()
    except RuntimeError as e:
        return _err(str(e))

    if not os.path.isfile(observation_csv):
        return _err(f"Observation file not found: {observation_csv}")

    # Extract "value" from each spec
    values = []
    for spec in param_specs:
        if "value" not in spec:
            return _err(
                f"param_spec is missing 'value' key: {spec}. "
                "Provide the value to test as spec['value']."
            )
        values.append(float(spec["value"]))

    t0 = time.time()

    # Resolve source directory for base case files
    project_dir = os.path.dirname(os.path.abspath(session.project_file))
    source_dir = os.path.abspath(base_case_dir) if base_case_dir else project_dir

    if not os.path.isdir(source_dir):
        return _err(f"base_case_dir not found: {source_dir}")

    # Create an isolated local run directory
    cleanup = run_dir is None
    if run_dir is None:
        parent = os.path.abspath(work_dir) if work_dir else os.getcwd()
        os.makedirs(parent, exist_ok=True)
        import tempfile as _tf
        run_dir = _tf.mkdtemp(prefix="calib_run_", dir=parent)

    # Log file for all pyHMT2D model output (print statements in model code)
    _model_log = os.path.join(
        os.path.abspath(work_dir) if work_dir else os.getcwd(), "pyHMT2D.log"
    )

    try:
        # Copy base case files and subdirectories into run_dir
        for fn in os.listdir(source_dir):
            src = os.path.join(source_dir, fn)
            dst = os.path.join(run_dir, fn)
            if os.path.isfile(src):
                shutil.copy2(src, dst)
            elif os.path.isdir(src):
                shutil.copytree(src, dst)

        # Re-open the copied project in the run_dir
        import pyHMT2D
        run_project_file = os.path.join(run_dir, os.path.basename(session.project_file))

        # Capture all pyHMT2D model print() output to log file + stderr.
        # This prevents stdout corruption of the MCP JSON protocol.
        with _capture_model_output(_model_log):
            if session.model_type == "SRH-2D":
                run_data = pyHMT2D.SRH_2D.SRH_2D_Data(run_project_file)
            else:
                # RAS_2D_Data expects the plan result HDF file (e.g. Muncie2D.p01.hdf),
                # not the .prj project file.  Derive its name from the session's known
                # hdf_filename (set when get_project_info opened the base case).
                if session.data is not None and hasattr(session.data, 'hdf_filename'):
                    base_hdf = os.path.basename(session.data.hdf_filename)
                else:
                    import glob as _glob
                    plan_hdfs = sorted(_glob.glob(os.path.join(run_dir, "*.p[0-9]*.hdf")))
                    if not plan_hdfs:
                        raise RuntimeError(
                            "No plan result HDF file found in run directory. "
                            "Open the project with get_project_info first."
                        )
                    base_hdf = os.path.basename(plan_hdfs[0])
                plan_hdf_abs = os.path.abspath(os.path.join(run_dir, base_hdf))
                if not os.path.isfile(plan_hdf_abs):
                    raise RuntimeError(
                        f"Plan result HDF not found in run directory: {plan_hdf_abs}"
                    )
                run_data = pyHMT2D.RAS_2D.RAS_2D_Data(plan_hdf_abs)

            # Apply parameter values
            _apply_params(run_data, session.model_type, param_specs, values)

            # Save modified inputs
            if session.model_type == "SRH-2D":
                ctrl = (
                    run_data.srhhydro_obj
                    if run_data.control_type == "SRHHydro"
                    else run_data.srhsif_obj
                )
                ctrl.save_as()

            # Run the model
            if session.model_type == "SRH-2D":
                exe = srh_path or getattr(session, "_srh_path", None)
                pre = srh_pre_path or getattr(session, "_srh_pre_path", exe)
                dll = extra_dll_path or getattr(session, "_extra_dll_path", "")
                if not exe:
                    raise RuntimeError("srh_path is required for SRH-2D calibration runs.")

                model = pyHMT2D.SRH_2D.SRH_2D_Model(version, pre, exe, dll, faceless=True)
                model.init_model()
                model.set_simulation_case(run_data)

                _orig_dir = os.getcwd()
                os.chdir(run_dir)
                try:
                    model.run_pre_model()
                    model.run_model(sleepTime=5.0, bShowProgress=False)
                finally:
                    os.chdir(_orig_dir)

                case_name = run_data.get_case_name()
                xmdf_file = os.path.join(run_dir, f"{case_name}_XMDFC.h5")
                if not os.path.isfile(xmdf_file):
                    import glob as _glob
                    candidates = _glob.glob(os.path.join(run_dir, "*.h5"))
                    xmdf_file = candidates[0] if candidates else xmdf_file

                run_data.readSRHXMDFFile(xmdf_file, bNodal=False)
                run_data.outputXMDFDataToVTK(bNodal=False, lastTimeStep=True, dir=run_dir)

            else:  # HEC-RAS
                import glob as _glob
                run_prj_files = _glob.glob(os.path.join(run_dir, "*.prj"))
                run_prj_file = run_prj_files[0] if run_prj_files else run_project_file
                model = pyHMT2D.RAS_2D.HEC_RAS_Model(version=hecras_version, faceless=True)
                model.init_model()
                model.open_project(run_prj_file)
                model.run_model()
                run_data = model.get_simulation_case()
                run_data.saveHEC_RAS2D_results_to_VTK(
                    lastTimeStep=True, fileNameBase="ras_result", dir=run_dir
                )
                model.close_project()
                model.exit_model()

        # Find produced VTK
        import glob as _glob
        vtk_files = sorted(_glob.glob(os.path.join(run_dir, "*.vtk")))
        if not vtk_files:
            return _err("Simulation produced no VTK output. Check solver logs.")
        vtk_file = vtk_files[-1]

        # Probe at observation points
        names, xs, ys, obs_vals = _parse_observation_csv(observation_csv)
        resolved_var = _resolve_result_variable(vtk_file, result_variable)
        sim_vals = _probe_wse_at_points(vtk_file, xs, ys, var_name=resolved_var)
        rmse = _rmse(sim_vals, obs_vals)

        per_point = [
            {
                "name": n,
                "x": x,
                "y": y,
                "simulated": round(s, 4),
                "measured": round(m, 4),
                "error": round(s - m, 4),
            }
            for n, x, y, s, m in zip(names, xs, ys, sim_vals, obs_vals)
        ]

        param_summary = {}
        for spec, val in zip(param_specs, values):
            key = f"{spec['type']}_{spec.get('material_id', spec.get('bc_id', '?'))}"
            param_summary[key] = val

        elapsed = round(time.time() - t0, 1)
        if cleanup and os.path.isdir(run_dir):
            shutil.rmtree(run_dir, ignore_errors=True)
        return _ok(
            {
                "rmse": round(rmse, 6),
                "per_point": per_point,
                "param_values": param_summary,
                "elapsed_seconds": elapsed,
                "run_dir": run_dir,
            },
            f"RMSE = {rmse:.4f} (elapsed {elapsed}s)",
        )

    except Exception as exc:
        # Leave the run directory intact on failure so logs can be inspected.
        return _err(
            f"Parameter evaluation failed: {exc}  |  run_dir: {run_dir}"
        )


# ── tool 22 ───────────────────────────────────────────────────────────────────

def run_calibration(
    param_specs: list,
    observation_csv: str,
    n_iterations: Optional[int] = None,
    method: str = "gp",
    result_variable: Optional[str] = None,
    base_case_dir: Optional[str] = None,
    work_dir: Optional[str] = None,
    output_csv: Optional[str] = None,
    srh_path: Optional[str] = None,
    srh_pre_path: Optional[str] = None,
    extra_dll_path: Optional[str] = None,
    version: str = "3.7",
    hecras_version: str = "6.6",
) -> dict:
    """Automatically calibrate model parameters to minimise RMSE against observations.

    Wraps an optimization loop around evaluate_parameters(). Each iteration
    copies the base case, applies a candidate parameter set, runs the model,
    and computes RMSE. The optimizer updates the next candidate based on the
    history of evaluations.

    Supports the same calibration workflow as the example scripts in
    examples/calibration/SRH-2D/ and examples/calibration/RAS-2D/.

    Parameters
    ----------
    param_specs : list of dict
        Calibration parameters. Each entry must have:
          {"type": "manning_n", "material_id": int, "material_name": str,
           "initial": float, "min": float, "max": float}
        The "value" key is ignored (the optimizer sets it).
    observation_csv : str
        Path to CSV with field observations (see load_observations()).
    n_iterations : int, optional
        Total number of model evaluations. If omitted, a sensible default is
        chosen automatically based on the number of calibration parameters:
        ``max(20, 10 * n_params + 5)`` for ``"gp"`` (Bayesian search needs
        enough points to build a useful surrogate), ``max(20, 5 * n_params)``
        for ``"nelder-mead"``.
    method : str
        Optimization algorithm:
          "gp"          — Gaussian Process Bayesian optimization (recommended,
                          from scikit-optimize gp_minimize). Global search,
                          works well with 2–5 parameters.
          "nelder-mead" — Nelder-Mead simplex (scipy.optimize). Local search,
                          faster but may find local minimum.
    result_variable : str
        VTK variable name to compare against observations. Default "Water_Elev_m".
    base_case_dir : str, optional
        Directory with the clean base-case files to copy from each iteration.
        Defaults to the directory of the open project file. A ``base_case/``
        subdirectory is the recommended layout.
    work_dir : str, optional
        Parent directory where per-iteration ``calib_run_XXXXXXXX`` subdirectories
        are created. Defaults to the current working directory.
    output_csv : str, optional
        Path to write the full evaluation history CSV. Defaults to
        'calibration_history.csv' in the project directory.
    srh_path, srh_pre_path, extra_dll_path : str
        SRH-2D solver paths (inherited from session if not provided).
    version : str
        SRH-2D version. Default "3.7".
    hecras_version : str
        HEC-RAS version. Default "6.6".

    Returns
    -------
    dict with "data" containing
        {"best_params": {name: value}, "best_rmse": float,
         "n_iterations": int, "history_csv": str, "elapsed_seconds": float}
    """
    session = get_session()
    try:
        session.require_open()
    except RuntimeError as e:
        return _err(str(e))

    if not os.path.isfile(observation_csv):
        return _err(f"Observation file not found: {observation_csv}")

    # Validate param_specs
    for spec in param_specs:
        for key in ("type", "initial", "min", "max"):
            if key not in spec:
                return _err(f"param_spec missing required key '{key}': {spec}")
        if spec["type"] == "manning_n" and "material_id" not in spec:
            return _err(f"Manning's n spec missing 'material_id': {spec}")

    # Output CSV path
    if output_csv is None:
        project_dir = os.path.dirname(os.path.abspath(session.project_file))
        output_csv = os.path.join(project_dir, "calibration_history.csv")

    n_params = len(param_specs)
    if n_iterations is None:
        if method == "gp":
            n_iterations = max(20, 10 * n_params + 5)
        else:  # nelder-mead
            n_iterations = max(20, 5 * n_params)

    # Progress log — written to work_dir so users can tail it in a terminal.
    _log_dir = os.path.abspath(work_dir) if work_dir else os.getcwd()
    os.makedirs(_log_dir, exist_ok=True)
    _log_path = os.path.join(_log_dir, "calib_progress.log")
    # Overwrite any existing log from a previous run.
    try:
        open(_log_path, "w").close()
    except OSError:
        pass

    def _log(msg: str) -> None:
        line = f"{time.strftime('%H:%M:%S')}  {msg}"
        print(line, file=_sys.stderr, flush=True)
        try:
            with open(_log_path, "a") as _fh:
                _fh.write(line + "\n")
        except OSError:
            pass

    _log(f"Calibration started: {n_iterations} iterations, method={method}, "
         f"{n_params} parameter(s)")

    bounds = [(float(spec["min"]), float(spec["max"])) for spec in param_specs]
    x0 = [float(spec["initial"]) for spec in param_specs]
    history = []
    t0 = time.time()

    # --- objective function ---
    def _objective(values):
        eval_num = len(history) + 1
        param_str = ", ".join(
            f"{s.get('material_name', s.get('type','p'))+'_n' if s['type']=='manning_n' else s['type']}={v:.4f}"
            for s, v in zip(param_specs, values)
        )
        _log(f"Evaluation {eval_num}/{n_iterations}  params: {param_str}")

        specs_with_vals = []
        for spec, val in zip(param_specs, values):
            s = dict(spec)
            s["value"] = float(val)
            specs_with_vals.append(s)

        result = evaluate_parameters(
            specs_with_vals,
            observation_csv,
            result_variable=result_variable,
            base_case_dir=base_case_dir,
            work_dir=work_dir,
            srh_path=srh_path,
            srh_pre_path=srh_pre_path,
            extra_dll_path=extra_dll_path,
            version=version,
            hecras_version=hecras_version,
        )
        if result["status"] == "ok":
            rmse = result["data"]["rmse"]
        else:
            rmse = 1e6
            _log(f"Evaluation {eval_num}/{n_iterations}  FAILED: {result['message']}")
        history.append({"values": list(values), "rmse": rmse})
        _log(f"Evaluation {eval_num}/{n_iterations}  RMSE = {rmse:.6f}")
        return rmse

    try:
        if method == "gp":
            from skopt import gp_minimize
            res = gp_minimize(
                _objective,
                dimensions=bounds,
                x0=x0,
                n_calls=n_iterations,
                random_state=42,
            )
            best_values = list(res.x)
            best_rmse = float(res.fun)

        elif method == "nelder-mead":
            from scipy.optimize import minimize

            res = minimize(
                _objective,
                x0,
                method="Nelder-Mead",
                bounds=bounds,
                options={"maxiter": n_iterations, "xatol": 1e-4, "fatol": 1e-6},
            )
            best_values = list(res.x)
            best_rmse = float(res.fun)

        else:
            return _err(
                f"Unknown method '{method}'. Supported: 'gp', 'nelder-mead'."
            )

        # Write history CSV
        param_names = [
            f"{s['type']}_{s.get('material_id', s.get('bc_id', i))}"
            for i, s in enumerate(param_specs)
        ]
        with open(output_csv, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(param_names + ["rmse"])
            for row in history:
                writer.writerow([round(v, 6) for v in row["values"]] + [round(row["rmse"], 6)])

        best_params = {
            name: round(val, 6)
            for name, val in zip(param_names, best_values)
        }
        elapsed = round(time.time() - t0, 1)
        best_str = ", ".join(f"{k}={v}" for k, v in best_params.items())
        _log(f"Calibration complete — best RMSE={best_rmse:.6f}  params: {best_str}  ({elapsed}s)")

        return _ok(
            {
                "best_params": best_params,
                "best_rmse": round(best_rmse, 6),
                "n_iterations_planned": n_iterations,
                "n_iterations_completed": len(history),
                "history_csv": output_csv,
                "elapsed_seconds": elapsed,
                "progress_log": _log_path,
            },
            f"Calibration complete. Best RMSE = {best_rmse:.4f} after "
            f"{len(history)} evaluations ({elapsed}s). "
            f"History saved to {output_csv}.",
        )

    except ImportError as e:
        return _err(
            f"Missing dependency for method='{method}': {e}. "
            "Install scikit-optimize (pip install scikit-optimize) for 'gp', "
            "or scipy (pip install scipy) for 'nelder-mead'."
        )
    except Exception as exc:
        return _err(f"Calibration failed: {exc}")
