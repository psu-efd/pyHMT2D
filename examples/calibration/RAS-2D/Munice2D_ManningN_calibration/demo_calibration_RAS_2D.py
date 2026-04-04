"""
RAS-2D calibration demo: calibrate Manning's n for channel (zone 2) and left_2 (zone 4)
using high water marks in HWMs.dat.

Uses only pyHMT2D basic functionality: HEC_RAS_Model, RAS_2D_Data, vtkHandler for
running HEC-RAS, converting results to VTK, and probing WSE at points. Optimization
is done with scikit-optimize's gp_minimize.

Run from this directory (the one containing base_case and HWMs.dat).
"""

import os
import subprocess
import sys
import shutil
import time
import numpy as np
from pathlib import Path

import pyHMT2D
from pyHMT2D.Misc import vtkHandler
import vtk

# Optional: scikit-optimize for Bayesian optimization
try:
    from skopt import gp_minimize
    from skopt.callbacks import CheckpointSaver
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
HEC_RAS_VERSION = "6.6"
PROJECT_FILE = "Muncie2D.prj"
TERRAIN_PATH = "Terrain/TerrainMuncie_composite.tif"
# Manning's n zone IDs (0-based in HEC-RAS): 1 = channel, 3 = left_2
MANNING_ZONE_IDS_0BASED = [1, 3]
MANNING_ZONE_NAMES = ["channel", "left_2"]
# Bounds for (n_channel, n_left_2)
BOUNDS = [(0.02, 0.08), (0.02, 0.08)]
N_CALLS = 10  # max objective evaluations. Adjust as needed (more calls = better calibration but longer runtime).


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
                data.append([parts[0], parts[1], parts[2]])  # x, y, WSE
    return np.array(data)


def run_hec_ras_with_mannings_n(run_dir, n_channel, n_left2):
    """
    Copy base_case to run_dir, set Manning's n for channel and left_2, run HEC-RAS,
    convert result to VTK. Returns (success, vtk_path or None).
    """
    run_dir = Path(run_dir)
    if run_dir.exists():
        # Retry rmtree — HEC-RAS may still hold file locks briefly after exit
        for attempt in range(5):
            try:
                shutil.rmtree(run_dir)
                break
            except PermissionError:
                if attempt < 4:
                    time.sleep(2)
                else:
                    # Kill any orphaned HEC-RAS processes locking the directory
                    print("Warning: calib_run still locked. Killing orphaned Ras.exe processes.")
                    subprocess.run(["taskkill", "/F", "/IM", "Ras.exe"],
                                   capture_output=True)
                    time.sleep(3)
                    shutil.rmtree(run_dir)  # final attempt
    shutil.copytree(BASE_CASE_DIR, run_dir)

    cwd = os.getcwd()
    model = None
    try:
        os.chdir(run_dir)

        model = pyHMT2D.RAS_2D.HEC_RAS_Model(version=HEC_RAS_VERSION, faceless=False)
        model.init_model()
        model.open_project(PROJECT_FILE)

        ras_data = model.get_simulation_case()
        # 0-based IDs for modify_ManningsN
        material_ids = [MANNING_ZONE_IDS_0BASED[0], MANNING_ZONE_IDS_0BASED[1]]
        ns = [n_channel, n_left2]
        ras_data.modify_ManningsN(material_ids, ns, MANNING_ZONE_NAMES)

        # Touch landcover so HEC-RAS picks up Manning change
        lc = ras_data.landcover_filename
        if lc and len(lc) > 0:
            lc_str = lc.decode("utf-8") if isinstance(lc, bytes) else lc
            lc_path = Path(run_dir) / lc_str if not Path(lc_str).is_absolute() else Path(lc_str)
            if lc_path.exists():
                lc_path.touch()

        model.save_project()
        ok = model.run_model()
        model.close_project()
        model.exit_model()
        model = None

        if not ok:
            os.chdir(cwd)
            return False, None

        # Convert to VTK using the new project-based API
        from pyHMT2D.Hydraulic_Models_Data.RAS_2D.HEC_RAS_Model import HEC_RAS_Project
        _project = HEC_RAS_Project(PROJECT_FILE)
        _plan = _project.get_plan("p01")
        ras_2d_data = _plan.load_results()
        vtk_list = ras_2d_data.saveHEC_RAS2D_results_to_VTK(lastTimeStep=True)
        vtk_path = Path(vtk_list[-1]).resolve()

        os.chdir(cwd)
        return True, vtk_path

    except Exception as e:
        print("Error in run_hec_ras_with_mannings_n:", e)
        if model is not None:
            try:
                model.close_project()
                model.exit_model()
            except Exception:
                pass
        os.chdir(cwd)
        return False, None


def probe_wse_at_points(vtk_path, points_xy):
    """
    Probe Water_Elev_ft at (x, y) points. points_xy: (N, 2) array.
    Returns (N,) array of WSE values.
    """
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
    x = [n_channel, n_left2].
    """
    global _evaluation_log
    n_channel, n_left2 = x[0], x[1]
    run_dir = SCRIPT_DIR / RUN_DIR_NAME

    success, vtk_path = run_hec_ras_with_mannings_n(run_dir, n_channel, n_left2)
    hwms = load_hwms(HWM_FILE)
    wse_measured = hwms[:, 2]

    if not success or vtk_path is None:
        _evaluation_log.append({
            "n_channel": n_channel, "n_left2": n_left2, "error": 1e10,
            "wse_sim": [np.nan] * len(wse_measured), "wse_meas": wse_measured.tolist(),
        })
        return 1e10  # large penalty on failure

    xy = hwms[:, :2]
    wse_sim = probe_wse_at_points(vtk_path, xy)
    if len(wse_sim) != len(wse_measured):
        _evaluation_log.append({
            "n_channel": n_channel, "n_left2": n_left2, "error": 1e10,
            "wse_sim": [np.nan] * len(wse_measured), "wse_meas": wse_measured.tolist(),
        })
        return 1e10

    rmse = np.sqrt(np.mean((wse_sim - wse_measured) ** 2))
    _evaluation_log.append({
        "n_channel": n_channel, "n_left2": n_left2, "error": rmse,
        "wse_sim": np.atleast_1d(wse_sim).tolist(), "wse_meas": wse_measured.tolist(),
    })

    # Clean run dir to save space
    if run_dir.exists():
        shutil.rmtree(run_dir, ignore_errors=True)

    return rmse


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


def calibrate_RAS_2D():
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
    print("Parameter bounds:", BOUNDS)
    print("Max evaluations:", N_CALLS)

    result = gp_minimize(
        objective,
        BOUNDS,
        n_calls=N_CALLS,
        random_state=42,
        verbose=True,
    )

    # Save calibration_results.csv (every evaluation), with WSE at each point
    n_points = len(hwms)
    header = "iteration,n_channel,n_left2,calibration_error"
    for p in range(n_points):
        header += f",wse_sim_point{p+1},wse_meas_point{p+1}"
    header += "\n"

    with open(RESULTS_CSV, "w") as f:
        f.write(header)
        for i, row in enumerate(_evaluation_log):
            line = f"{i+1},{row['n_channel']},{row['n_left2']},{row['error']}"
            for p in range(n_points):
                sim = row["wse_sim"][p] if p < len(row["wse_sim"]) else np.nan
                meas = row["wse_meas"][p] if p < len(row["wse_meas"]) else np.nan
                line += f",{sim},{meas}"
            line += "\n"
            f.write(line)

    print("Best parameters: n_channel =", result.x[0], ", n_left2 =", result.x[1])
    print("Best RMSE =", result.fun)
    print("Results written to", RESULTS_CSV.name)


    # Copy base_case to calibrated_case and run with best parameters
    if CALIBRATED_CASE_DIR.exists():
        shutil.rmtree(CALIBRATED_CASE_DIR)
    shutil.copytree(BASE_CASE_DIR, CALIBRATED_CASE_DIR)

    print("Running calibrated case and saving to", CALIBRATED_CASE_DIR.name, "...")
    success, vtk_path = run_hec_ras_with_mannings_n(
        CALIBRATED_CASE_DIR,
        result.x[0],
        result.x[1],
    )
    if success:
        print("Calibrated case VTK:", vtk_path)
    else:
        print("Calibrated case run had errors.")

    print("Done.")


if __name__ == "__main__":
    """Run the calibration."""
    calibrate_RAS_2D()

    # Plot calibration trajectory on objective contour (read from CSV)
    plot_calibration_trajectory(
        RESULTS_CSV,
        BOUNDS,
        save_path=SCRIPT_DIR / "calibration_trajectory.png",
    )


