"""
Note: This does not work well. The L-BFGS-B optimizer does not converge. Maybe the gradient is too noisy.


SRH-2D calibration demo: calibrate Manning's n for selected material zones
using high water marks in HWMs.dat.

Uses pyHMT2D (SRH_2D_Model, SRH_2D_Data, vtkHandler) to run SRH-2D, convert
results to VTK, and probe WSE at HWM points. Optimization is done with
scipy.optimize.minimize (L-BFGS-B for bounded parameters).

Run from this directory (the one containing base_case and HWMs.dat).
HWMs.dat format: x, y, WSE (measured), one point per line; lines starting with # are skipped.
"""

import os
import sys
import shutil
import numpy as np
from pathlib import Path

import pyHMT2D
from pyHMT2D.Misc import vtkHandler
import vtk
from scipy.optimize import minimize

# --- Configuration (modify for your case) ---
SCRIPT_DIR = Path(__file__).resolve().parent
BASE_CASE_DIR = SCRIPT_DIR / "base_case"
HWM_FILE = SCRIPT_DIR / "HWMs.dat"
RUN_DIR_NAME = "calib_run"
CALIBRATED_CASE_DIR = SCRIPT_DIR / "calibrated_case"
RESULTS_CSV = SCRIPT_DIR / "calibration_results.csv"

# SRH-2D control file (inside base_case / run dir)
SRH_CONTROL_FILE = "Muncie.srhhydro"

# SRH-2D executable paths (edit for your system)
VERSION = "3.6.5"
SRH_PRE_PATH = r"C:\Program Files\SMS 13.4 64-bit\python\Lib\site-packages\srh2d_exe\SRH_Pre_Console.exe"
SRH_PATH = r"C:\Program Files\SMS 13.4 64-bit\python\Lib\site-packages\srh2d_exe\SRH-2D_Console.exe"
EXTRA_DLL_PATH = r"C:\Program Files\SMS 13.4 64-bit\python\Lib\site-packages\srh2d_exe"

# Manning's n: material IDs to calibrate (ID from srhhydro ManningsN lines)
# Base case has 0..7 with values 0.03, 0.03, 0.03, 0.04, 0.045, 0.05, 0.055, 0.06
MANNING_MATERIAL_IDS = [2, 4]
MANNING_MATERIAL_NAMES = ["channel", "left_2"]
# Bounds for (n_zone2, n_zone3)
BOUNDS = [(0.02, 0.08), (0.02, 0.08)]
# Initial guess 
X0 = [0.03, 0.04]
# Optimizer options
MAXITER = 20
FTOL = 1e-6
# Step size for numerical gradient (finite differences). Default 1e-8 is too small
# for Manning's n (~0.03): model output doesn't change, gradient is zero, optimizer stops.
# Use a fraction of typical parameter scale so RMSE actually changes when n is perturbed.
EPS_GRAD = 1e-3  # absolute perturbation ~0.001 for gradient; increase if optimizer still doesn't move


def load_hwms(filepath):
    """Load high water marks from file: x, y, WSE (skip comment lines)."""
    data = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [float(x) for x in line.split(",")]
            if len(parts) >= 3:
                data.append([parts[0], parts[1], parts[2]])
    return np.array(data)


def run_srh_2d_with_mannings_n(run_dir, ns):
    """
    Copy base_case to run_dir, set Manning's n for the calibrated materials,
    run SRH-2D Pre and model, convert result to VTK.
    ns: list of Manning n values (same length as MANNING_MATERIAL_IDS).
    Returns (success, vtk_path or None).
    """
    run_dir = Path(run_dir)
    if run_dir.exists():
        shutil.rmtree(run_dir)
    shutil.copytree(BASE_CASE_DIR, run_dir)

    cwd = os.getcwd()
    try:
        os.chdir(run_dir)

        my_srh_2d_model = pyHMT2D.SRH_2D.SRH_2D_Model(
            VERSION, SRH_PRE_PATH, SRH_PATH, EXTRA_DLL_PATH, faceless=False
        )
        my_srh_2d_model.init_model()
        my_srh_2d_model.open_project(SRH_CONTROL_FILE)

        my_srh_2d_data = my_srh_2d_model.get_simulation_case()
        my_srh_2d_data.modify_ManningsNs(
            MANNING_MATERIAL_IDS, ns, MANNING_MATERIAL_NAMES
        )
        my_srh_2d_data.srhhydro_obj.save_as()

        ok = my_srh_2d_model.run_pre_model()
        if not ok:
            my_srh_2d_model.close_project()
            my_srh_2d_model.exit_model()
            os.chdir(cwd)
            return False, None

        ok = my_srh_2d_model.run_model()
        my_srh_2d_model.close_project()
        my_srh_2d_model.exit_model()

        if not ok:
            os.chdir(cwd)
            return False, None

        # Convert XMDF to VTK
        my_srh_2d_data = pyHMT2D.SRH_2D.SRH_2D_Data(SRH_CONTROL_FILE)
        bNodal = False
        xmdf_file = my_srh_2d_data.get_case_name() + "_XMDFC.h5"
        my_srh_2d_data.readSRHXMDFFile(xmdf_file, bNodal)
        vtk_list = my_srh_2d_data.outputXMDFDataToVTK(bNodal, lastTimeStep=True, dir="")
        vtk_path = Path(vtk_list[-1]).resolve()

        os.chdir(cwd)
        return True, vtk_path

    except Exception as e:
        print("Error in run_srh_2d_with_mannings_n:", e)
        os.chdir(cwd)
        return False, None


def probe_wse_at_points(vtk_path, points_xy):
    """Probe Water_Elev_ft at (x, y) points. points_xy: (N, 2) array. Returns (N,) WSE."""
    vtk_handler_obj = vtkHandler()
    reader = vtk_handler_obj.readVTK_UnstructuredGrid(str(vtk_path))

    pt_vtk = vtk.vtkPoints()
    pt_vtk.SetNumberOfPoints(len(points_xy))
    for i, (x, y) in enumerate(points_xy):
        pt_vtk.InsertPoint(i, x, y, 0.0)

    points, WSE, elev = vtk_handler_obj.probeUnstructuredGridVTKOnPoints(
        pt_vtk, reader, "Water_Elev_ft"
    )
    return np.atleast_1d(WSE)


# Log every objective evaluation for calibration_results.csv
_evaluation_log = []


def objective(x):
    """
    Objective for minimize: RMSE between simulated and measured WSE at HWM points.
    x = [n_channel, n_left_2] (or however many MANNING_MATERIAL_IDS).
    """
    global _evaluation_log
    ns = list(x)
    run_dir = SCRIPT_DIR / RUN_DIR_NAME

    success, vtk_path = run_srh_2d_with_mannings_n(run_dir, ns)
    hwms = load_hwms(HWM_FILE)
    wse_measured = hwms[:, 2]

    if not success or vtk_path is None:
        row = {**{f"n{i}": ns[i] for i in range(len(ns))}, "error": 1e10,
               "wse_sim": [np.nan] * len(wse_measured), "wse_meas": wse_measured.tolist()}
        _evaluation_log.append(row)
        n_evals = len(_evaluation_log)
        param_str = ", ".join(
            f"{MANNING_MATERIAL_NAMES[i]}={ns[i]:.4f}" for i in range(len(ns))
        )
        print(f"  Evaluation {n_evals}: {param_str}  RMSE = failed")
        return 1e10

    xy = hwms[:, :2]
    wse_sim = probe_wse_at_points(vtk_path, xy)
    if len(wse_sim) != len(wse_measured):
        row = {**{f"n{i}": ns[i] for i in range(len(ns))}, "error": 1e10,
               "wse_sim": [np.nan] * len(wse_measured), "wse_meas": wse_measured.tolist()}
        _evaluation_log.append(row)
        n_evals = len(_evaluation_log)
        param_str = ", ".join(
            f"{MANNING_MATERIAL_NAMES[i]}={ns[i]:.4f}" for i in range(len(ns))
        )
        print(f"  Evaluation {n_evals}: {param_str}  RMSE = failed")
        return 1e10

    rmse = np.sqrt(np.mean((wse_sim - wse_measured) ** 2))
    row = {**{f"n{i}": ns[i] for i in range(len(ns))}, "error": rmse,
           "wse_sim": np.atleast_1d(wse_sim).tolist(), "wse_meas": wse_measured.tolist()}
    _evaluation_log.append(row)
    # Print progress (callback is not invoked by L-BFGS-B in scipy)
    n_evals = len(_evaluation_log)
    param_str = ", ".join(
        f"{MANNING_MATERIAL_NAMES[i]}={ns[i]:.4f}" for i in range(len(ns))
    )
    print(f"  Evaluation {n_evals}: {param_str}  RMSE = {rmse:.4f}")

    if run_dir.exists():
        shutil.rmtree(run_dir, ignore_errors=True)

    return rmse


def main():
    global _evaluation_log
    _evaluation_log = []

    if not BASE_CASE_DIR.is_dir():
        print("base_case directory not found. Run from the calibration example directory.")
        sys.exit(1)
    if not HWM_FILE.is_file():
        print("HWMs.dat not found.")
        sys.exit(1)

    hwms = load_hwms(HWM_FILE)
    print("Loaded", len(hwms), "high water marks from", HWM_FILE.name)
    print("Calibrating Manning's n for material IDs:", MANNING_MATERIAL_IDS)
    print("Bounds:", BOUNDS)
    print("Initial guess:", X0)
    print("Method: L-BFGS-B, maxiter =", MAXITER, ", eps (gradient step) =", EPS_GRAD)
    print("Progress (each line = one objective evaluation, including gradient finite-diff):")

    # L-BFGS-B uses finite differences for the gradient when jac is not provided.
    # eps: step size for those finite differences. Default 1e-8 is too small for
    # Manning's n (~0.03), so the objective doesn't change and the optimizer stops.
    # Note: L-BFGS-B in scipy does not invoke callback; progress is printed from objective().
    result = minimize(
        objective,
        X0,
        method="L-BFGS-B",
        bounds=BOUNDS,
        options=dict(maxiter=MAXITER, ftol=FTOL, eps=EPS_GRAD),
    )

    # Build header and save calibration_results.csv
    n_params = len(MANNING_MATERIAL_IDS)
    n_points = len(hwms)
    param_names = [f"n{MANNING_MATERIAL_NAMES[i]}" for i in range(n_params)]
    header = "iteration," + ",".join(param_names) + ",calibration_error"
    for p in range(n_points):
        header += f",wse_sim_point{p+1},wse_meas_point{p+1}"
    header += "\n"

    with open(RESULTS_CSV, "w") as f:
        f.write(header)
        for i, row in enumerate(_evaluation_log):
            line = f"{i+1}," + ",".join(str(row[f"n{j}"]) for j in range(n_params)) + f",{row['error']}"
            for p in range(n_points):
                sim = row["wse_sim"][p] if p < len(row["wse_sim"]) else np.nan
                meas = row["wse_meas"][p] if p < len(row["wse_meas"]) else np.nan
                line += f",{sim},{meas}"
            line += "\n"
            f.write(line)

    print("Optimization converged:", result.success)
    print("Best parameters:", dict(zip(MANNING_MATERIAL_NAMES, result.x)))
    print("Best RMSE =", result.fun)
    print("Results written to", RESULTS_CSV.name)

    # Run calibrated case and save to calibrated_case
    if CALIBRATED_CASE_DIR.exists():
        shutil.rmtree(CALIBRATED_CASE_DIR)
    shutil.copytree(BASE_CASE_DIR, CALIBRATED_CASE_DIR)

    print("Running calibrated case and saving to", CALIBRATED_CASE_DIR.name, "...")
    success, vtk_path = run_srh_2d_with_mannings_n(CALIBRATED_CASE_DIR, list(result.x))
    if success:
        print("Calibrated case VTK:", vtk_path)
    else:
        print("Calibrated case run had errors.")

    print("Done.")


if __name__ == "__main__":
    main()
