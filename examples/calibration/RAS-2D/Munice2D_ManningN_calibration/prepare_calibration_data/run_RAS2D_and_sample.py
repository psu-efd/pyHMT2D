"""
Run HEC-RAS 2D, convert the results to VTK, and sample the results on VTK at the high water marks.

High Water Marks: x, y
  413253.87, 1800957.98
  411429.26, 1803300.59

Uses pyHMT2D to control HEC-RAS 2D (HEC_RAS_Model), convert results to VTK (RAS_2D_Data),
and probe WSE at HWM points (vtkHandler). Writes sampled WSE to HWMs.dat.

Run from this directory (where Muncie2D.prj and Terrain/ are).
"""

import os
import sys
from pathlib import Path

import numpy as np
import pyHMT2D
from pyHMT2D.Misc import vtkHandler
import vtk

# --- Configuration (modify for your case) ---
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_FILE = "Muncie2D.prj"
TERRAIN_PATH = "Terrain/TerrainMuncie_composite.tif"
PLAN_HDF = "Muncie2D.p01.hdf"
OUTPUT_SAMPLED = SCRIPT_DIR / "HWMs.dat"

# High water marks: x, y
HWM_POINTS = np.array([
    [413253.87, 1800957.98],
    [411429.26, 1803300.59],
])

HEC_RAS_VERSION = "6.6"


def run_hec_ras_2d_and_export_vtk():
    """
    Run HEC-RAS 2D and export last time step to VTK.
    Returns (success, vtk_file_path or None).
    """
    model = pyHMT2D.RAS_2D.HEC_RAS_Model(version=HEC_RAS_VERSION, faceless=False)
    model.init_model()
    model.open_project(PROJECT_FILE)

    ok = model.run_model()
    model.close_project()
    model.exit_model()

    if not ok:
        print("HEC-RAS 2D run failed.")
        return False, None

    # Convert to VTK
    ras_2d_data = pyHMT2D.RAS_2D.RAS_2D_Data(PLAN_HDF, TERRAIN_PATH)
    vtk_list = ras_2d_data.saveHEC_RAS2D_results_to_VTK(lastTimeStep=True)
    vtk_path = Path(vtk_list[-1]).resolve()
    return True, vtk_path


def probe_wse_at_points(vtk_path, points_xy, var_name="Water_Elev_ft"):
    """
    Probe WSE at (x, y) points. points_xy: (N, 2) array.
    Returns (N,) array of WSE values.
    """
    vtk_handler_obj = vtkHandler()
    reader = vtk_handler_obj.readVTK_UnstructuredGrid(str(vtk_path))

    pt_vtk = vtk.vtkPoints()
    pt_vtk.SetNumberOfPoints(len(points_xy))
    for i, (x, y) in enumerate(points_xy):
        pt_vtk.InsertPoint(i, x, y, 0.0)

    points, WSE, elev = vtk_handler_obj.probeUnstructuredGridVTKOnPoints(
        pt_vtk, reader, var_name
    )
    return np.atleast_1d(WSE)


def main():
    os.chdir(SCRIPT_DIR)

    print("Running HEC-RAS 2D and exporting to VTK...")
    ok, vtk_path = run_hec_ras_2d_and_export_vtk()
    if not ok or vtk_path is None:
        sys.exit(1)

    print("Probing WSE at high water marks...")
    points_xy = HWM_POINTS
    simulated_wse = probe_wse_at_points(vtk_path, points_xy)

    with open(OUTPUT_SAMPLED, "w") as f:
        f.write("# x, y, simulated_WSE\n")
        for i in range(len(points_xy)):
            f.write(
                f"{points_xy[i, 0]}, {points_xy[i, 1]}, {simulated_wse[i]}\n"
            )

    print("Sampled results written to", OUTPUT_SAMPLED)
    print("  Point 1: simulated WSE =", simulated_wse[0])
    print("  Point 2: simulated WSE =", simulated_wse[1])


if __name__ == "__main__":
    main()
    print("Done.")
