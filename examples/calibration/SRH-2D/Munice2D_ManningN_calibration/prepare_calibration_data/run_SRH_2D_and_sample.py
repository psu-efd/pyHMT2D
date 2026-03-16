"""
Run SRH-2D, convert the results to VTK, and sample the results on VTK at the high water marks.

High Water Marks: x, y
  413253.87, 1800957.98
  411429.26, 1803300.59

Uses pyHMT2D to control SRH-2D (SRH_2D_Model), convert XMDF to VTK (SRH_2D_Data),
and probe WSE at HWM points (vtkHandler). Writes sampled WSE to HWMs_sampled.dat.

Run from this directory (where Muncie.srhhydro, Muncie.srhgeom, Muncie.srhmat are).
"""

import os
import sys
from pathlib import Path

import numpy as np
import pyHMT2D
from pyHMT2D.Misc import vtkHandler
import vtk

# --- Paths (modify for your SRH-2D installation) ---
SCRIPT_DIR = Path(__file__).resolve().parent
SRH_CONTROL_FILE = "Muncie.srhhydro"
OUTPUT_SAMPLED = SCRIPT_DIR / "HWMs.dat"
# High water marks: x, y
HWM_POINTS = np.array([
    [413253.87, 1800957.98],
    [411429.26, 1803300.59],
])

# SRH-2D executable paths (edit for your system)
VERSION = "3.6.5"
SRH_PRE_PATH = r"C:\Program Files\SMS 13.4 64-bit\python\Lib\site-packages\srh2d_exe\SRH_Pre_Console.exe"
SRH_PATH = r"C:\Program Files\SMS 13.4 64-bit\python\Lib\site-packages\srh2d_exe\SRH-2D_Console.exe"
EXTRA_DLL_PATH = r"C:\Program Files\SMS 13.4 64-bit\python\Lib\site-packages\srh2d_exe"


def run_srh_2d_and_export_vtk():
    """
    Run SRH-2D and export last time step to VTK.
    Returns (success, vtk_file_path or None).
    """
    my_srh_2d_model = pyHMT2D.SRH_2D.SRH_2D_Model(
        VERSION, SRH_PRE_PATH, SRH_PATH, EXTRA_DLL_PATH, faceless=False
    )
    my_srh_2d_model.init_model()

    print("Hydraulic model:", my_srh_2d_model.getName(), my_srh_2d_model.getVersion())

    my_srh_2d_model.open_project(SRH_CONTROL_FILE)
    ok = my_srh_2d_model.run_pre_model()
    if not ok:
        print("SRH-2D Pre failed.")
        my_srh_2d_model.close_project()
        my_srh_2d_model.exit_model()
        return False, None

    ok = my_srh_2d_model.run_model()
    my_srh_2d_model.close_project()
    my_srh_2d_model.exit_model()

    if not ok:
        print("SRH-2D run failed.")
        return False, None

    # Convert XMDF result to VTK
    my_srh_2d_data = pyHMT2D.SRH_2D.SRH_2D_Data(SRH_CONTROL_FILE)
    bNodal = False
    xmdf_file = my_srh_2d_data.get_case_name() + "_XMDFC.h5"
    my_srh_2d_data.readSRHXMDFFile(xmdf_file, bNodal)
    vtk_list = my_srh_2d_data.outputXMDFDataToVTK(bNodal, lastTimeStep=True, dir="")
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

    print("Running SRH-2D and exporting to VTK...")
    ok, vtk_path = run_srh_2d_and_export_vtk()
    if not ok or vtk_path is None:
        sys.exit(1)

    print("Probing WSE at high water marks...")
    points_xy = HWM_POINTS[:, :2]
    simulated_wse = probe_wse_at_points(vtk_path, points_xy)

    # Write x, y, measured_WSE, simulated_WSE
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
