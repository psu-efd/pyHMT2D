# -*- coding: utf-8 -*-
"""
Group 8 — Format conversion tools.

These tools are stateless: they take explicit file paths as arguments and
produce output files. No open project or session is required.
"""

from __future__ import annotations

import os


def _ok(data, message: str) -> dict:
    return {"status": "ok", "data": data, "message": message}


def _err(message: str) -> dict:
    return {"status": "error", "data": None, "message": message}


# ── tool 27 ───────────────────────────────────────────────────────────────────

def ras_to_srh(
    ras_hdf_file: str,
    terrain_tif_file: str,
    srh_case_name: str,
) -> dict:
    """Convert a HEC-RAS 2D mesh and Manning's n values to SRH-2D format.

    Reads the HEC-RAS 2D result HDF file and terrain GeoTIFF, then writes
    ``<srh_case_name>.srhgeom`` and ``<srh_case_name>.srhmat`` in the current
    working directory.  The ``.srhhydro`` boundary-condition file must be
    created separately.

    Parameters
    ----------
    ras_hdf_file : str
        HEC-RAS 2D result HDF file (e.g. ``Muncie2D.p01.hdf``).
    terrain_tif_file : str
        Terrain file in GeoTIFF format (e.g. ``Terrain/terrain.tif``).
    srh_case_name : str
        Base name for the output SRH-2D files (e.g. ``srh_Muncie`` →
        produces ``srh_Muncie.srhgeom`` and ``srh_Muncie.srhmat``).

    Returns
    -------
    dict with "data" containing
        {"srhgeom_file": str, "srhmat_file": str}
    """
    if not os.path.isfile(ras_hdf_file):
        return _err(f"HEC-RAS HDF file not found: {ras_hdf_file}")
    if not os.path.isfile(terrain_tif_file):
        return _err(f"Terrain GeoTIFF not found: {terrain_tif_file}")

    try:
        import pyHMT2D
        converter = pyHMT2D.Misc.RAS_to_SRH_Converter(
            ras_hdf_file, terrain_tif_file, srh_case_name
        )
        converter.convert_to_SRH()

        srhgeom = f"{srh_case_name}.srhgeom"
        srhmat = f"{srh_case_name}.srhmat"
        return _ok(
            {"srhgeom_file": srhgeom, "srhmat_file": srhmat},
            f"Converted HEC-RAS to SRH-2D: {srhgeom}, {srhmat}",
        )
    except Exception as exc:
        return _err(f"Conversion failed: {exc}")


# ── tool 28 ───────────────────────────────────────────────────────────────────

def srh_to_vtk(
    srhhydro_file: str,
    output_file: str,
    mesh_only: bool = False,
    nodal: bool = False,
    all_timesteps: bool = False,
) -> dict:
    """Convert SRH-2D results to VTK, or export only the mesh.

    Parameters
    ----------
    srhhydro_file : str
        SRH-2D case file in ``.srhhydro`` format.
    output_file : str
        When ``mesh_only=False``: path to the SRH-2D HDF result file
        (``.h5``); VTK output file(s) are written alongside it.
        When ``mesh_only=True``: path for the output ``.vtk`` mesh file.
    mesh_only : bool
        If True, export only the computational mesh (no result data).
        Default False.
    nodal : bool
        If True, use nodal (vertex) data instead of cell-centred data.
        Ignored when ``mesh_only=True``. Default False.
    all_timesteps : bool
        If True, export all time steps to separate VTK files; otherwise
        only the last time step is exported. Ignored when
        ``mesh_only=True``. Default False.

    Returns
    -------
    dict with "data" containing {"vtk_files": list[str]}
    """
    if not os.path.isfile(srhhydro_file):
        return _err(f"SRH-2D hydro file not found: {srhhydro_file}")
    if not mesh_only and not os.path.isfile(output_file):
        return _err(f"SRH-2D HDF result file not found: {output_file}")

    try:
        import pyHMT2D
        srh_data = pyHMT2D.SRH_2D.SRH_2D_Data(srhhydro_file)

        if mesh_only:
            srh_data.srhgeom_obj.output_2d_mesh_to_vtk(output_file, bFlat=False)
            return _ok(
                {"vtk_files": [output_file]},
                f"Exported SRH-2D mesh to {output_file}",
            )
        else:
            srh_data.readSRHXMDFFile(output_file, bNodal=nodal)
            srh_data.outputXMDFDataToVTK(
                bNodal=nodal,
                lastTimeStep=not all_timesteps,
                dir='',
            )
            # Collect generated VTK files in cwd
            vtk_files = sorted(
                f for f in os.listdir('.') if f.endswith('.vtk')
            )
            return _ok(
                {"vtk_files": vtk_files},
                f"Exported SRH-2D results to VTK: {len(vtk_files)} file(s)",
            )
    except Exception as exc:
        return _err(f"Export failed: {exc}")
