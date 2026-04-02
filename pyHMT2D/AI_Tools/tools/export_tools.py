# -*- coding: utf-8 -*-
"""
Group 5 — Export tools.

Convert simulation results to VTK format for visualization in ParaView,
and export the mesh geometry.
"""

from __future__ import annotations

import os
import tempfile
from typing import Optional

from pyHMT2D.AI_Tools.state import get_session


def _ok(data, message: str) -> dict:
    return {"status": "ok", "data": data, "message": message}


def _err(message: str) -> dict:
    return {"status": "error", "data": None, "message": message}


# ── tool 18 ───────────────────────────────────────────────────────────────────

def export_to_vtk(
    timestep: int = -1,
    output_dir: Optional[str] = None,
    nodal: bool = False,
) -> dict:
    """Export simulation results to VTK format for visualization in ParaView.

    For SRH-2D: converts the loaded XMDF result to one or more .vtk files.
    For HEC-RAS: converts the loaded HDF result to .vtk files.
    Requires results to be loaded first (call read_results()).

    Parameters
    ----------
    timestep : int
        Time step index to export. Use -1 to export the last time step only.
        Use -2 to export all time steps (may produce many files).
        Default -1.
    output_dir : str, optional
        Directory to write VTK files to. Defaults to a subdirectory named
        'vtk_output' next to the project file.
    nodal : bool
        (SRH-2D only) If True, export nodal (vertex) data instead of
        cell-centred data. Default False (cell-centred).

    Returns
    -------
    dict with "data" containing
        {"vtk_files": [str, ...], "output_dir": str}
    """
    session = get_session()
    try:
        session.require_results()
    except RuntimeError as e:
        return _err(str(e))

    # Resolve output directory
    if output_dir is None:
        project_dir = (
            os.path.dirname(session.project_file) if session.project_file else "."
        )
        output_dir = os.path.join(project_dir, "vtk_output")
    os.makedirs(output_dir, exist_ok=True)

    try:
        import glob

        if session.model_type == "SRH-2D":
            last = timestep == -1
            all_steps = timestep == -2
            ts = timestep if not (last or all_steps) else -1

            session.data.outputXMDFDataToVTK(
                bNodal=nodal,
                timeStep=ts,
                lastTimeStep=last,
                dir=output_dir,
            )

        else:  # HEC-RAS
            last = timestep == -1
            ts = timestep if not last else -1
            session.data.saveHEC_RAS2D_results_to_VTK(
                timeStep=ts,
                lastTimeStep=last,
                fileNameBase="ras_result",
                dir=output_dir,
                bFlat=False,
            )

        vtk_files = sorted(glob.glob(os.path.join(output_dir, "*.vtk")))

        # Update session VTK cache to the last file produced
        if vtk_files:
            session.vtk_file = vtk_files[-1]

        return _ok(
            {"vtk_files": vtk_files, "output_dir": output_dir},
            f"Exported {len(vtk_files)} VTK file(s) to {output_dir}.",
        )

    except Exception as exc:
        return _err(f"VTK export failed: {exc}")


# ── tool 19 ───────────────────────────────────────────────────────────────────

def export_mesh_to_vtk(
    output_path: Optional[str] = None,
    flat: bool = False,
) -> dict:
    """Export the computational mesh (without results) to a VTK file.

    Useful for inspecting mesh quality, cell counts, and material zone layout
    before running a simulation. Requires a project to be open.

    Parameters
    ----------
    output_path : str, optional
        Full path for the output VTK file. Defaults to 'mesh.vtk' in a
        'vtk_output' subdirectory next to the project file.
    flat : bool
        If True, set all z-coordinates to 0 (useful for 2D visualisation).
        Default False (preserves bed elevation).

    Returns
    -------
    dict with "data" containing {"vtk_file": str}
    """
    session = get_session()
    try:
        session.require_open()
    except RuntimeError as e:
        return _err(str(e))

    if session.model_type != "SRH-2D":
        return _err(
            "export_mesh_to_vtk() currently supports SRH-2D only. "
            "For HEC-RAS mesh, use export_to_vtk() after loading results."
        )

    if output_path is None:
        project_dir = (
            os.path.dirname(session.project_file) if session.project_file else "."
        )
        vtk_dir = os.path.join(project_dir, "vtk_output")
        os.makedirs(vtk_dir, exist_ok=True)
        output_path = os.path.join(vtk_dir, "mesh.vtk")

    try:
        session.data.output_2d_mesh_to_vtk(
            meshVTKFileName=os.path.basename(output_path),
            bFlat=flat,
            dir=os.path.dirname(output_path),
        )

        return _ok(
            {"vtk_file": output_path},
            f"Mesh exported to {output_path}.",
        )

    except Exception as exc:
        return _err(f"Mesh export failed: {exc}")
