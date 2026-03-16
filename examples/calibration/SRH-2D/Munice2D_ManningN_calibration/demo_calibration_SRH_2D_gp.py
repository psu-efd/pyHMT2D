"""
SRH-2D calibration demo: calibrate Manning's n for selected material zones
using high water marks in HWMs.dat.

Uses pyHMT2D (SRH_2D_Model, SRH_2D_Data, vtkHandler) to run SRH-2D, convert
results to VTK, and probe WSE at HWM points. Optimization is done with
scikit-optimize's gp_minimize (Gaussian Process Bayesian optimization).

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

try:
    from skopt import gp_minimize
    HAS_SKOPT = True
except ImportError:
    HAS_SKOPT = False

try:
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.rcParams["font.family"] = "serif"
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    
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
VERSION = "3.7.1"
SRH_PRE_PATH = r"C:\Program Files\SMS 13.4 64-bit\python\Lib\site-packages\srh2d_exe\SRH_Pre_Console.exe"
SRH_PATH = r"C:\Program Files\SMS 13.4 64-bit\python\Lib\site-packages\srh2d_exe\SRH-2D_Console.exe"
EXTRA_DLL_PATH = r"C:\Program Files\SMS 13.4 64-bit\python\Lib\site-packages\srh2d_exe"

# Manning's n: material IDs to calibrate (ID from srhhydro ManningsN lines)
MANNING_MATERIAL_IDS = [2, 4]
MANNING_MATERIAL_NAMES = ["channel", "left_2"]
# Bounds for each parameter (same order as MANNING_MATERIAL_IDS)
BOUNDS = [(0.02, 0.08), (0.02, 0.08)]
N_CALLS = 50  # max objective evaluations (gp_minimize)


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

        ok = my_srh_2d_model.run_model(sleepTime=5.0)
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
    Objective for gp_minimize: RMSE between simulated and measured WSE at HWM points.
    x = [n_param1, n_param2, ...] (same length as MANNING_MATERIAL_IDS).
    """
    global _evaluation_log
    ns = list(x)
    run_dir = SCRIPT_DIR / RUN_DIR_NAME

    success, vtk_path = run_srh_2d_with_mannings_n(run_dir, ns)
    hwms = load_hwms(HWM_FILE)
    wse_measured = hwms[:, 2]

    if not success or vtk_path is None:
        row = {
            **{f"n{i}": ns[i] for i in range(len(ns))},
            "error": 1e10,
            "wse_sim": [np.nan] * len(wse_measured),
            "wse_meas": wse_measured.tolist(),
        }
        _evaluation_log.append(row)
        n_evals = len(_evaluation_log)
        param_str = ", ".join(
            f"{MANNING_MATERIAL_NAMES[i]}={ns[i]:.4f}" for i in range(len(ns))
        )
        print(f"  Iteration {n_evals}/{N_CALLS}: {param_str}  RMSE = failed")
        if run_dir.exists():
            shutil.rmtree(run_dir, ignore_errors=True)
        return 1e10

    xy = hwms[:, :2]
    wse_sim = probe_wse_at_points(vtk_path, xy)
    if len(wse_sim) != len(wse_measured):
        row = {
            **{f"n{i}": ns[i] for i in range(len(ns))},
            "error": 1e10,
            "wse_sim": [np.nan] * len(wse_measured),
            "wse_meas": wse_measured.tolist(),
        }
        _evaluation_log.append(row)
        n_evals = len(_evaluation_log)
        param_str = ", ".join(
            f"{MANNING_MATERIAL_NAMES[i]}={ns[i]:.4f}" for i in range(len(ns))
        )
        print(f"  Iteration {n_evals}/{N_CALLS}: {param_str}  RMSE = failed")
        if run_dir.exists():
            shutil.rmtree(run_dir, ignore_errors=True)
        return 1e10

    rmse = np.sqrt(np.mean((wse_sim - wse_measured) ** 2))
    row = {
        **{f"n{i}": ns[i] for i in range(len(ns))},
        "error": rmse,
        "wse_sim": np.atleast_1d(wse_sim).tolist(),
        "wse_meas": wse_measured.tolist(),
    }
    _evaluation_log.append(row)
    n_evals = len(_evaluation_log)
    param_str = ", ".join(
        f"{MANNING_MATERIAL_NAMES[i]}={ns[i]:.4f}" for i in range(len(ns))
    )
    print(f"  Iteration {n_evals}/{N_CALLS}: {param_str}  RMSE = {rmse:.4f}")

    if run_dir.exists():
        shutil.rmtree(run_dir, ignore_errors=True)

    return rmse


def calibrate_SRH_2D():
    """Calibrate the Manning's n for the channel and left_2 in the Muncie 2D case."""

    global _evaluation_log
    _evaluation_log = []

    if not BASE_CASE_DIR.is_dir():
        print("base_case directory not found. Run from the calibration example directory.")
        sys.exit(1)
    if not HWM_FILE.is_file():
        print("HWMs.dat not found.")
        sys.exit(1)
    if not HAS_SKOPT:
        print("Install scikit-optimize: pip install scikit-optimize")
        sys.exit(1)

    hwms = load_hwms(HWM_FILE)
    print("Loaded", len(hwms), "high water marks from", HWM_FILE.name)
    print("Calibrating Manning's n for material IDs:", MANNING_MATERIAL_IDS)
    print("Parameter bounds:", BOUNDS)
    print("Max evaluations (n_calls):", N_CALLS)
    print("Method: gp_minimize (Gaussian Process Bayesian optimization)")
    print("Progress (each line = one objective evaluation):")

    result = gp_minimize(
        objective,
        BOUNDS,
        n_calls=N_CALLS,
        random_state=42,
        verbose=True,
    )

    # Save calibration_results.csv
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

def plot_calibration_trajectory(csv_path, bounds, save_path=None):
    """
    Plot calibration trajectory on top of objective function contours (n_channel, n_left2).
    Reads evaluation history from calibration_results.csv; trajectory in iteration order;
    best parameter point (min calibration_error) is highlighted and labeled.
    Publication-ready formatting.
    """
    if not HAS_MATPLOTLIB:
        print("matplotlib not available; skipping calibration trajectory plot.")
        return

    csv_path = Path(csv_path)
    if not csv_path.is_file():
        print("Calibration results CSV not found:", csv_path)
        return

    # Read from CSV: iteration, n_channel, n_left2, calibration_error, ...
    data = np.loadtxt(csv_path, delimiter=",", skiprows=1, usecols=(1, 2, 3))
    n_ch = data[:, 0]
    n_lt = data[:, 1]
    err = data[:, 2]

    # Best point = row with minimum calibration_error
    best_idx = np.argmin(err)
    x_best = (float(n_ch[best_idx]), float(n_lt[best_idx]))
    fun_best = float(err[best_idx])

    # Use finite/successful evaluations for contour (cap penalty so contour is sensible)
    penalty = 1e10
    finite = err < penalty * 0.1
    if not np.any(finite):
        # Fallback: scatter only
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.scatter(n_ch, n_lt, c="gray", s=20, alpha=0.7, label="Evaluations")
        ax.scatter([x_best[0]], [x_best[1]], marker="*", s=400, color="gold", edgecolor="black", zorder=5, label="Optimized")
        ax.set_xlim(bounds[0])
        ax.set_ylim(bounds[1])
        ax.set_xlabel(r"$n_{\mathrm{channel}}$", fontsize=12)
        ax.set_ylabel(r"$n_{\mathrm{left\_2}}$", fontsize=12)
        ax.legend(loc="lower right", fontsize=10)
        ax.set_aspect("equal")
        plt.tight_layout()
        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()
        return

    err_plot = np.where(finite, err, np.nan)
    err_finite = err[finite]
    n_ch_f = n_ch[finite]
    n_lt_f = n_lt[finite]
    vmax = np.percentile(err_finite, 95) if len(err_finite) > 2 else np.max(err_finite)
    vmin = np.min(err_finite)

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    # Filled contours from scattered data (only successful points for smoothness)
    try:
        levels = np.linspace(vmin, vmax, 14)
        tc = ax.tricontourf(n_ch_f, n_lt_f, err_finite, levels=levels, cmap="viridis_r", alpha=0.85)
        ax.tricontour(n_ch_f, n_lt_f, err_finite, levels=levels, colors="k", linewidths=0.3, alpha=0.4)
    except Exception:
        # Fallback if triangulation fails (e.g. too few points)
        tc = ax.scatter(n_ch_f, n_lt_f, c=err_finite, s=80, cmap="viridis_r", alpha=0.9)

    cbar = fig.colorbar(tc, ax=ax, shrink=0.8, aspect=25)
    cbar.set_label(r"RMSE (ft)", fontsize=11)
    cbar.ax.tick_params(labelsize=9)

    # Tried parameter points (no connecting lines)
    ax.scatter(n_ch, n_lt, c=np.arange(len(n_ch)), cmap="plasma", s=32, edgecolor="white", linewidths=0.5, zorder=4, label="Evaluations")
    # Best point
    ax.scatter(
        [x_best[0]], [x_best[1]],
        marker="*", s=450, color="gold", edgecolor="black", linewidths=1.2, zorder=5, label="Optimized",
    )
    label = (
        r"Optimized: $n_{\mathrm{ch}}=%.3f$, $n_{\mathrm{left2}}=%.3f$" + "\n" + r"RMSE = %.4f ft"
    ) % (x_best[0], x_best[1], fun_best)
    ax.annotate(
        label,
        xy=(x_best[0], x_best[1]),
        xytext=(10, 18),
        textcoords="offset points",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="wheat", alpha=0.95, edgecolor="gray"),
        arrowprops=dict(arrowstyle="->", color="black", lw=1.0),
        zorder=6,
    )

    ax.set_xlim(bounds[0])
    ax.set_ylim(bounds[1])
    ax.set_xlabel(r"$n_{\mathrm{channel}}$", fontsize=12)
    ax.set_ylabel(r"$n_{\mathrm{left\_2}}$", fontsize=12)
    ax.tick_params(axis="both", labelsize=10)
    ax.set_aspect("equal")
    ax.legend(loc="lower right", fontsize=10)
    ax.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print("Calibration trajectory plot saved to", save_path.name)
    plt.close()

    print("Done.")

if __name__ == "__main__":
    """Run the calibration."""
    calibrate_SRH_2D()
    
    # Plot calibration trajectory on objective contour (read from CSV)
    plot_calibration_trajectory(
        RESULTS_CSV,
        BOUNDS,
        save_path=SCRIPT_DIR / "calibration_trajectory.png",
    )


