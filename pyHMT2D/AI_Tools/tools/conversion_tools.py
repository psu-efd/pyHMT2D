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
    srh_case_name: str,
    terrain_tif_file: str = "",
) -> dict:
    """Convert a HEC-RAS 2D mesh and Manning's n values to SRH-2D format.

    Derives the HEC-RAS project (.prj) file automatically from the HDF path,
    then writes ``<srh_case_name>.srhgeom`` and ``<srh_case_name>.srhmat`` in
    the current working directory.  The ``.srhhydro`` boundary-condition file
    must be created separately.

    Parameters
    ----------
    ras_hdf_file : str
        HEC-RAS 2D result HDF file (e.g. ``Muncie2D.p01.hdf``). Used to
        derive the project file and plan ID automatically.
    srh_case_name : str
        Base name for the output SRH-2D files (e.g. ``srh_Muncie`` →
        produces ``srh_Muncie.srhgeom`` and ``srh_Muncie.srhmat``).
    terrain_tif_file : str, optional
        Kept for backwards compatibility; ignored (terrain is auto-detected
        from the project geometry HDF and .rasmap file).

    Returns
    -------
    dict with "data" containing
        {"srhgeom_file": str, "srhmat_file": str}
    """
    if not os.path.isfile(ras_hdf_file):
        return _err(f"HEC-RAS HDF file not found: {ras_hdf_file}")

    try:
        import re as _re
        import glob as _glob

        # Derive .prj file from the HDF path (e.g. Muncie2D.p01.hdf -> Muncie2D.prj)
        hdf_abs = os.path.abspath(ras_hdf_file)
        hdf_dir = os.path.dirname(hdf_abs)
        stem = os.path.basename(hdf_abs)
        for ext in (".hdf", ".hdf5"):
            if stem.lower().endswith(ext):
                stem = stem[: -len(ext)]
        # Extract plan ID (e.g. "p01") and strip it from the stem
        plan_match = _re.search(r'\.([pP]\d+)$', stem)
        plan_id = plan_match.group(1).lower() if plan_match else None
        stem_no_plan = _re.sub(r'\.[pP]\d+$', '', stem)
        candidate_prj = os.path.join(hdf_dir, stem_no_plan + ".prj")
        if not os.path.isfile(candidate_prj):
            prj_files = _glob.glob(os.path.join(hdf_dir, "*.prj"))
            if not prj_files:
                return _err(
                    f"Could not find a .prj project file next to: {ras_hdf_file}. "
                    "Pass the .prj file path as ras_hdf_file or place the .prj alongside the HDF."
                )
            candidate_prj = prj_files[0]

        if plan_id is None:
            return _err(
                f"Could not determine plan ID from HDF filename: {os.path.basename(ras_hdf_file)}. "
                "Expected a name like 'CaseName.p01.hdf'."
            )

        import pyHMT2D
        converter = pyHMT2D.Misc.RAS_to_SRH_Converter(
            candidate_prj, plan_id, srh_case_name
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
