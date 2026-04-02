# -*- coding: utf-8 -*-
"""
Group 1 — Discovery tools.

These are the first tools an AI agent calls to understand what model files
exist and what parameters / areas are in the currently open project.
All tools return a dict with keys:
    "status"  : "ok" | "error"
    "data"    : the payload (list / dict / scalar)
    "message" : human-readable summary
"""

from __future__ import annotations

import json
import os
import glob as _glob
from typing import Optional

import pyHMT2D
from pyHMT2D.AI_Tools.state import get_session


# ── helpers ───────────────────────────────────────────────────────────────────

def _ok(data, message: str) -> dict:
    return {"status": "ok", "data": data, "message": message}


def _err(message: str) -> dict:
    return {"status": "error", "data": None, "message": message}


def _find_prj_for_hdf(hdf_path: str) -> Optional[str]:
    """Derive the HEC-RAS .prj file path from a plan result .hdf path.

    Example: ``Muncie2D.p01.hdf`` → ``Muncie2D.prj`` (same directory).
    Strips the plan/geometry shortID suffix (e.g. ``.p01``, ``.g01``) from
    the stem before appending ``.prj``.  Falls back to any ``.prj`` found in
    the same directory if the derived name does not exist.
    """
    import re

    d = os.path.dirname(os.path.abspath(hdf_path))
    stem = os.path.basename(hdf_path)
    for ext in (".hdf", ".hdf5"):
        if stem.lower().endswith(ext):
            stem = stem[: -len(ext)]
    # strip plan/geometry shortID like .p01, .g01, .u01 …
    stem = re.sub(r"\.[pgu]\d+$", "", stem, flags=re.IGNORECASE)
    candidate = os.path.join(d, stem + ".prj")
    if os.path.isfile(candidate):
        return candidate
    # fallback: any .prj in the directory
    prj_files = _glob.glob(os.path.join(d, "*.prj"))
    return prj_files[0] if prj_files else None


def _detect_model_type(file_path: str) -> Optional[str]:
    """Infer model type from file extension and (for .prj) content."""
    lower = file_path.lower()
    if lower.endswith(".srhhydro") or lower.endswith("_sif.dat") or lower.endswith(".sif"):
        return "SRH-2D"
    if lower.endswith(".hdf") or lower.endswith(".hdf5"):
        return "HEC-RAS"
    if lower.endswith(".prj"):
        # A HEC-RAS .prj file contains "HEC-RAS Version" on the first few lines.
        try:
            with open(file_path, "r", errors="ignore") as fh:
                head = fh.read(512)
            if "HEC-RAS" in head or "Proj Title" in head:
                return "HEC-RAS"
        except OSError:
            pass
    return None


# ── tool 1 ────────────────────────────────────────────────────────────────────

def list_model_files(directory: str) -> dict:
    """List hydraulic model input and result files in a directory.

    Searches for SRH-2D control files (.srhhydro, _SIF.dat), HEC-RAS project
    files (.prj), and HEC-RAS result files (.hdf). Returns each file labelled
    with its model type so the AI agent can pick the right one to open.

    Parameters
    ----------
    directory : str
        Path to the directory to search (non-recursive).

    Returns
    -------
    dict with "data" containing a list of
        {"file": str, "model_type": str, "case_name": str}
    """
    if not os.path.isdir(directory):
        return _err(f"Directory not found: {directory}")

    import re as _re

    # HEC-RAS plan result HDF files follow the naming convention <case>.p<nn>.hdf
    # (e.g. Muncie2D.p01.hdf).  Geometry (.g01.hdf), terrain, Manning's n, and
    # other auxiliary HDF files must be excluded — they are not project files.
    _ras_plan_hdf_re = _re.compile(r"\.[pP]\d+\.hdf5?$")

    srh2d_patterns = ["*.srhhydro", "*_SIF.dat", "*.sif"]
    # For HEC-RAS list .prj files and plan-result HDF files only
    ras_patterns   = ["*.prj", "*.hdf", "*.hdf5"]

    found = []

    # SRH-2D files
    for ext in srh2d_patterns:
        for fp in _glob.glob(os.path.join(directory, ext)):
            if _detect_model_type(fp) != "SRH-2D":
                continue
            case_name = os.path.splitext(os.path.basename(fp))[0]
            if case_name.upper().endswith("_SIF"):
                case_name = case_name[:-4]
            found.append({"file": fp, "model_type": "SRH-2D", "case_name": case_name})

    # HEC-RAS files: .prj (project) and *.p<nn>.hdf (plan results) only
    for ext in ras_patterns:
        for fp in _glob.glob(os.path.join(directory, ext)):
            name = os.path.basename(fp)
            lower = name.lower()
            if lower.endswith(".prj"):
                if _detect_model_type(fp) != "HEC-RAS":
                    continue
                case_name = os.path.splitext(name)[0]
            elif _ras_plan_hdf_re.search(name):
                # Plan result HDF — keep; strip .p01.hdf suffix for case name
                case_name = _ras_plan_hdf_re.sub("", name)
            else:
                continue  # geometry, terrain, Manning's n, or other auxiliary HDF
            found.append({"file": fp, "model_type": "HEC-RAS", "case_name": case_name})

    if not found:
        return _ok([], f"No hydraulic model files found in {directory}.")

    found.sort(key=lambda x: x["file"])
    return _ok(found, f"Found {len(found)} model file(s) in {directory}.")


# ── tool 2 ────────────────────────────────────────────────────────────────────

def get_project_info(project_file: str) -> dict:
    """Open a hydraulic model project and return its metadata.

    Auto-detects the model type from the file extension. Stores the opened
    model and data objects in the session so all subsequent tools can use them
    without re-specifying the model type.

    For SRH-2D: accepts a .srhhydro or _SIF.dat control file.
    For HEC-RAS: accepts a .prj project file or a .hdf result file.

    Parameters
    ----------
    project_file : str
        Path to the SRH-2D control file or HEC-RAS project / HDF file.

    Returns
    -------
    dict with "data" containing project metadata (case name, units, time info,
    materials count, boundary count, etc.)
    """
    if not os.path.isfile(project_file):
        return _err(f"File not found: {project_file}")

    model_type = _detect_model_type(project_file)
    if model_type is None:
        return _err(
            f"Cannot determine model type from file: {project_file}. "
            "Expected .srhhydro, _SIF.dat, .prj, or .hdf."
        )

    session = get_session()

    try:
        if model_type == "SRH-2D":
            data = pyHMT2D.SRH_2D.SRH_2D_Data(project_file)
            case_name = data.get_case_name()
            output_fmt, output_unit = data.get_output_format_unit()
            n_dict = data.get_ManningN_dict()
            n_zones = len(n_dict) if n_dict else 0

            # Boundary conditions
            if data.control_type == "SRHHydro":
                bc_dict = data.srhhydro_obj.get_BC()
            else:
                bc_dict = data.srhsif_obj.srhsif_content.get("BC", {})
            n_bcs = len(bc_dict) if bc_dict else 0

            # Mesh stats
            n_cells = data.srhgeom_obj.numOfElements if data.srhgeom_obj else None
            n_nodes = data.srhgeom_obj.numOfNodes if data.srhgeom_obj else None

            info = {
                "model_type": "SRH-2D",
                "project_file": project_file,
                "case_name": case_name,
                "control_type": data.control_type,
                "output_format": output_fmt,
                "output_unit": output_unit,
                "n_material_zones": n_zones,
                "n_boundary_conditions": n_bcs,
                "n_cells": n_cells,
                "n_nodes": n_nodes,
            }

            session.model_type = "SRH-2D"
            session.project_file = project_file
            session.data = data
            session.model = None  # model object created on run_simulation

        else:  # HEC-RAS
            # Accept either .prj (preferred) or a plan result .hdf.
            # Always resolve to the .prj file so COM can open the project.
            if project_file.lower().endswith(".prj"):
                prj_file = os.path.abspath(project_file)
            else:
                prj_file = _find_prj_for_hdf(project_file)
                if prj_file is None:
                    return _err(
                        f"Could not find a .prj project file next to: {project_file}. "
                        "Pass the .prj file directly."
                    )

            if not os.path.isfile(prj_file):
                return _err(f"HEC-RAS project file not found: {prj_file}")

            # HEC-RAS version from session config (set by load_config), default "6.6"
            ras_version = getattr(session, "_hecras_version", None) or "6.6"

            # Open via COM — requires HEC-RAS to be installed
            try:
                model = pyHMT2D.RAS_2D.HEC_RAS_Model(version=ras_version, faceless=True)
            except Exception as exc:
                return _err(
                    f"Failed to create HEC_RAS_Model (is HEC-RAS {ras_version} installed?): {exc}"
                )

            try:
                current_plan_name, current_plan_file = model.open_project(prj_file)
            except Exception as exc:
                return _err(f"Failed to open HEC-RAS project via COM: {exc}")

            proj = model._project

            # Build plan summary list
            plans_info = []
            for p in proj.plans:
                plan_hdf = (p.plan_file + ".hdf") if p.plan_file else None
                plans_info.append({
                    "name": p.plan_name,
                    "plan_file": p.plan_file,
                    "geom_file": p.geom_file,
                    "flow_file": p.flow_file,
                    "steady": p.steady,
                    "result_hdf": plan_hdf,
                    "result_hdf_exists": os.path.isfile(plan_hdf) if plan_hdf else False,
                })

            # Optionally load RAS_2D_Data for the current plan's result HDF
            # (only available after the simulation has been run at least once)
            current_plan_hdf = (current_plan_file + ".hdf") if current_plan_file else None
            data = None
            n_zones = 0
            n_bcs = None
            area_names = []
            cell_counts = []
            hecras_version_str = ""
            units = ""
            start_time = None
            end_time = None
            plan_short_id = ""

            if current_plan_hdf and os.path.isfile(current_plan_hdf):
                hdf_abs = os.path.abspath(current_plan_hdf)
                hdf_dir = os.path.dirname(hdf_abs)
                _orig_dir = os.getcwd()
                os.chdir(hdf_dir)
                try:
                    data = pyHMT2D.RAS_2D.RAS_2D_Data(os.path.basename(hdf_abs))
                finally:
                    os.chdir(_orig_dir)
                # Pin hdf_filename to absolute path so h5py works from any CWD
                data.hdf_filename = hdf_abs

                hecras_version_str = data.version
                units = data.units
                start_time = data.start_time
                end_time = data.end_time
                # Decode bytes area names and convert numpy ints to plain Python ints
                raw_areas = data.TwoDAreaNames or []
                area_names = [
                    a.decode("utf-8", errors="replace") if isinstance(a, bytes) else str(a)
                    for a in raw_areas
                ]
                raw_counts = data.TwoDAreaCellCounts or []
                cell_counts = [int(c) for c in raw_counts]
                n_zones = len(data.ManningNZones) if data.ManningNZones else 0
                plan_short_id = data.plan_shortID
                try:
                    n_bcs = len(data.TwoDAreaBoundaryNamesTypes)
                except Exception:
                    n_bcs = None

            n_cells = int(sum(cell_counts)) if cell_counts else None

            info = {
                "model_type": "HEC-RAS",
                "project_file": prj_file,
                "project_title": proj.title,
                "current_plan_name": current_plan_name,
                "current_plan_file": current_plan_file,
                "plan_shortID": plan_short_id,
                "hecras_version": hecras_version_str,
                "units": units,
                "simulation_start": str(start_time) if start_time else None,
                "simulation_end": str(end_time) if end_time else None,
                "2d_flow_areas": area_names,
                "cell_counts_per_area": cell_counts,
                "n_cells": n_cells,
                "n_material_zones": n_zones,
                "n_boundary_conditions": n_bcs,
                "plans": plans_info,
                "terrain_hdf_files": proj.terrain_hdf_file_list,
                "terrain_tiff_files": proj.terrain_tiff_file_list,
                "result_hdf_available": data is not None,
            }

            session.model_type = "HEC-RAS"
            session.project_file = prj_file
            session._hecras_prj_file = prj_file
            session.data = data   # None if no result HDF yet
            session.model = model  # keep COM connection alive

        return _ok(info, f"Opened {model_type} project: {project_file}")

    except Exception as exc:
        return _err(f"Failed to open project: {exc}")


# ── tool 3 ────────────────────────────────────────────────────────────────────

def get_materials() -> dict:
    """Return all material zones and their current Manning's n values.

    Requires a project to be open (call get_project_info first).
    Use the returned material IDs and names when calling set_manning_n().

    Returns
    -------
    dict with "data" containing a list of
        {"id": int, "name": str, "manning_n": float}
    """
    session = get_session()
    try:
        session.require_open()
    except RuntimeError as e:
        return _err(str(e))

    try:
        if session.model_type == "SRH-2D":
            n_dict = session.data.get_ManningN_dict()
            # n_dict is {material_id (int): manning_n (float), ...}
            # Material names are in srhmat if available
            # Build mat_names: {int_id: name} from SRH_2D_SRHMat.matNameList.
            # matNameList stores string keys ("1", "2", …, "-1" for default).
            mat_names = {}
            if session.data.srhmat_obj is not None:
                raw = getattr(session.data.srhmat_obj, "matNameList", None) or {}
                for k, v in raw.items():
                    try:
                        mat_names[int(k)] = str(v).strip('"').strip("'")
                    except (ValueError, TypeError):
                        pass

            materials = []
            for mat_id, n_val in sorted(n_dict.items()):
                # Try int key first, then string key as fallback
                name = (
                    mat_names.get(mat_id)
                    or mat_names.get(str(mat_id))
                    or f"Material_{mat_id}"
                )
                materials.append({
                    "id": mat_id,
                    "name": name,
                    "manning_n": float(n_val),
                })

        else:  # HEC-RAS
            # ManningNZones: {zone_id: [name_bytes, n_value], ...}
            if not session.data.ManningNZones:
                session.data.build2DManningNZones()
            materials = []
            for zone_id, zone_info in sorted(session.data.ManningNZones.items()):
                name = zone_info[0]
                if isinstance(name, bytes):
                    name = name.decode("utf-8", errors="replace")
                materials.append({
                    "id": zone_id,
                    "name": name,
                    "manning_n": float(zone_info[1]),
                })

        return _ok(
            materials,
            f"Found {len(materials)} material zone(s) in {session.model_type} project.",
        )

    except Exception as exc:
        return _err(f"Failed to get materials: {exc}")


# ── tool 4 ────────────────────────────────────────────────────────────────────

def get_boundary_conditions() -> dict:
    """Return all boundary conditions with their IDs, types, and current values.

    Requires a project to be open (call get_project_info first).
    Use the returned IDs when calling set_inlet_flow() or set_exit_wse().

    Returns
    -------
    dict with "data" containing a list of
        {"id": int/str, "name": str, "type": str, "value": float or None}
    """
    session = get_session()
    try:
        session.require_open()
    except RuntimeError as e:
        return _err(str(e))

    try:
        bcs = []

        if session.model_type == "SRH-2D":
            ctrl = (
                session.data.srhhydro_obj
                if session.data.control_type == "SRHHydro"
                else session.data.srhsif_obj
            )

            if session.data.control_type == "SRHHydro":
                bc_dict = ctrl.get_BC()
                iq_dict = ctrl.get_InletQ()
                exh_dict = ctrl.get_ExitH()
            else:
                bc_dict = ctrl.srhsif_content.get("BC", {})
                iq_dict = ctrl.srhsif_content.get("IQParams", {})
                exh_dict = ctrl.srhsif_content.get("EWSParamsC", {})

            _, output_unit = session.data.get_output_format_unit()
            is_english = output_unit.upper() in ("EN", "ENGLISH")
            q_units = "ft3/s" if is_english else "m3/s"
            h_units = "ft"   if is_english else "m"

            for bc_id, bc_type in sorted(bc_dict.items()):
                value = None
                if bc_type == "INLET-Q" and bc_id in iq_dict:
                    value = float(iq_dict[bc_id].get("discharge", 0))
                elif bc_type == "EXIT-H" and bc_id in exh_dict:
                    entry = exh_dict[bc_id]
                    value = float(entry) if not isinstance(entry, dict) else float(entry.get("wse", 0))
                bcs.append({
                    "id": bc_id,
                    "name": f"BC_{bc_id}",
                    "type": bc_type,
                    "value": value,
                    "value_units": q_units if bc_type == "INLET-Q" else (h_units if bc_type == "EXIT-H" else None),
                })

        else:  # HEC-RAS
            bc_info = session.data.get2DAreaBoundaryNamesTypes()
            for i, entry in enumerate(bc_info):
                name = entry[0].decode("utf-8", errors="replace") if isinstance(entry[0], bytes) else str(entry[0])
                area = entry[1].decode("utf-8", errors="replace") if isinstance(entry[1], bytes) else str(entry[1])
                btype = entry[2].decode("utf-8", errors="replace") if isinstance(entry[2], bytes) else str(entry[2])
                length = float(entry[3]) if entry[3] is not None else None
                bcs.append({
                    "id": i,
                    "name": name,
                    "2d_flow_area": area,
                    "type": btype,
                    "length_m": length,
                    "value": None,  # HEC-RAS BC values live in plan files
                })

        return _ok(bcs, f"Found {len(bcs)} boundary condition(s).")

    except Exception as exc:
        return _err(f"Failed to get boundary conditions: {exc}")


# ── tool 5 ────────────────────────────────────────────────────────────────────

def get_result_variables(result_file: str) -> dict:
    """Return available result variables and time steps from a result file.

    For SRH-2D: accepts an XMDF/HDF5 result file (.h5, .xmdf).
    For HEC-RAS: the result file is the .hdf plan output file.

    This tool does NOT require a project to be open first — it can be called
    directly to inspect a result file. However, if a project is already open
    its model type is used; otherwise the model type is inferred from the file.

    Parameters
    ----------
    result_file : str
        Path to the SRH-2D XMDF file or HEC-RAS HDF result file.

    Returns
    -------
    dict with "data" containing
        {"variables": [...], "n_time_steps": int, "time_values": [...]}
    """
    if not os.path.isfile(result_file):
        return _err(f"Result file not found: {result_file}")

    session = get_session()

    try:
        import h5py

        model_type = session.model_type  # may be None

        # Infer if session not open
        if model_type is None:
            lower = result_file.lower()
            if lower.endswith(".h5") or lower.endswith(".xmdf"):
                model_type = "SRH-2D"
            elif lower.endswith(".hdf") or lower.endswith(".hdf5"):
                model_type = "HEC-RAS"
            else:
                return _err(
                    "Cannot determine model type from result file extension. "
                    "Open a project first with get_project_info()."
                )

        variables = []
        time_values = []

        if model_type == "SRH-2D":
            # SRH-2D XMDF structure: variables are top-level groups.
            # Each group has "Times" and "Values" datasets.
            # e.g.  /Water_Elev_m/Times, /Water_Elev_m/Values
            _SKIP = {"File Type", "File Version"}
            with h5py.File(result_file, "r") as hf:
                for key in hf.keys():
                    if key in _SKIP:
                        continue
                    grp = hf[key]
                    if isinstance(grp, h5py.Group) and "Values" in grp:
                        variables.append(key)
                        # Grab time values from the first variable found
                        if not time_values and "Times" in grp:
                            time_values = grp["Times"][:].tolist()

        else:  # HEC-RAS
            with h5py.File(result_file, "r") as hf:
                try:
                    results_grp = hf["Results"]["Unsteady"]["Output"]["Output Blocks"][
                        "Base Output"
                    ]["Unsteady Time Series"]["2D Flow Areas"]
                    area_name = list(results_grp.keys())[0]
                    area_grp = results_grp[area_name]
                    variables = [
                        k for k in area_grp.keys()
                        if isinstance(area_grp[k], h5py.Dataset)
                    ]
                    # Times
                    try:
                        t_key = hf["Results"]["Unsteady"]["Output"]["Output Blocks"][
                            "Base Output"
                        ]["Unsteady Time Series"]["Time Date Stamp"]
                        time_values = [
                            t.decode("utf-8", errors="replace") if isinstance(t, bytes) else str(t)
                            for t in t_key[:]
                        ]
                    except KeyError:
                        pass
                except KeyError:
                    pass

        # Deduplicate and sort
        variables = sorted(set(variables))
        n_steps = len(time_values)

        # For HEC-RAS: also compute the VTK variable names that query tools
        # (get_result_statistics, get_value_at_point) actually use.  The VTK
        # names are derived from the project units and differ from the HDF
        # names returned above (e.g. "Water Surface" → "Water_Elev_ft").
        vtk_variables = None
        if model_type == "HEC-RAS":
            session = get_session()
            units = getattr(session.data, "units", None) if session.data else None
            if units == "Feet":
                len_sfx, vel_sfx = "_ft", "_ft_p_s"
            elif units in ("Meter", "Metres", "Meters"):
                len_sfx, vel_sfx = "_m", "_m_p_s"
            else:
                len_sfx, vel_sfx = None, None
            if len_sfx:
                vtk_variables = sorted([
                    f"Water_Elev{len_sfx}",
                    f"Water_Depth{len_sfx}",
                    f"Velocity_cell{vel_sfx}",
                    "ManningN",
                    f"Velocity_node{vel_sfx}",
                    f"Bed_Elev{len_sfx}",
                ])

        info = {
            "model_type": model_type,
            "result_file": result_file,
            "variables": variables,
            "n_time_steps": n_steps,
            "time_values": time_values[:20],  # cap at 20 for readability
            "time_values_truncated": n_steps > 20,
        }
        if vtk_variables is not None:
            info["query_variable_names"] = vtk_variables
            info["note"] = (
                "Use 'query_variable_names' (not 'variables') when calling "
                "get_result_statistics(), get_value_at_point(), get_flood_extent(), "
                "get_cross_section_profile(), or export_to_vtk()."
            )

        msg = f"Found {len(variables)} variable(s) over {n_steps} time step(s)."
        if vtk_variables:
            msg += f" Use query_variable_names for result queries."
        return _ok(info, msg)

    except Exception as exc:
        return _err(f"Failed to read result variables: {exc}")


# ── helper tool ───────────────────────────────────────────────────────────────

def build_param_specs(specs: list) -> dict:
    """Validate and build a param_specs list for calibration or Monte Carlo tools.

    Translates a user-friendly parameter description into the exact dict format
    expected by ``run_calibration()``, ``evaluate_parameters()``,
    ``generate_mc_samples()``, and ``run_monte_carlo()``.

    The AI agent should call this tool after ``get_materials()`` /
    ``get_boundary_conditions()`` so it can resolve material names to IDs and
    confirm every entry is valid before starting a long run.

    Each entry in ``specs`` must have:

    - ``"type"`` : ``"manning_n"`` | ``"inlet_q"`` | ``"exit_h"``
    - For ``"manning_n"``: either ``"material_id"`` (int) **or**
      ``"material_name"`` (str, matched case-insensitively to open project).
    - For ``"inlet_q"`` / ``"exit_h"``: ``"bc_id"`` (int).
    - ``"min"``, ``"max"`` : parameter bounds (required for all types).
    - ``"initial"`` : starting value for calibration (optional; defaults to
      midpoint of ``[min, max]``).
    - ``"distribution"`` : ``"truncated_normal"`` | ``"uniform"`` (Monte Carlo
      only; default ``"truncated_normal"``).
    - ``"mean"``, ``"std"`` : distribution parameters for ``truncated_normal``
      (default to midpoint and ``(max-min)/6``).

    Parameters
    ----------
    specs : list of dict
        One entry per calibration/MC parameter (see above).

    Returns
    -------
    dict with ``"data"`` containing the validated and completed param_specs list,
    ready to pass directly to calibration or Monte Carlo tools.
    """
    session = get_session()
    try:
        session.require_open()
    except RuntimeError as e:
        return _err(str(e))

    # Build name→id lookup from open project materials
    mat_name_to_id: dict = {}
    mat_id_to_name: dict = {}
    try:
        if session.model_type == "SRH-2D":
            raw = getattr(session.data.srhmat_obj, "matNameList", None) or {}
            for k, v in raw.items():
                try:
                    mid = int(k)
                except (ValueError, TypeError):
                    continue
                clean = str(v).strip('"').strip("'")
                mat_name_to_id[clean.lower()] = mid
                mat_id_to_name[mid] = clean
        else:  # HEC-RAS: ManningNZones = {id: [name, n_value]}
            if session.data is not None:
                zones = getattr(session.data, "ManningNZones", None) or {}
                for mid, zone_info in zones.items():
                    name = zone_info[0]
                    if isinstance(name, bytes):
                        name = name.decode("utf-8", errors="replace")
                    name = str(name).strip()
                    mat_id_to_name[int(mid)] = name
                    mat_name_to_id[name.lower()] = int(mid)
    except Exception:
        pass  # name lookup is best-effort; ID-based specs still work

    errors = []
    result_specs = []

    for i, spec in enumerate(specs):
        s = dict(spec)
        label = f"spec[{i}]"

        ptype = s.get("type")
        if ptype not in ("manning_n", "inlet_q", "exit_h"):
            errors.append(f"{label}: 'type' must be 'manning_n', 'inlet_q', or 'exit_h' (got {ptype!r}).")
            continue

        # Resolve material_id from name if needed
        if ptype == "manning_n":
            if "material_id" not in s:
                name = s.get("material_name", "")
                matched = mat_name_to_id.get(name.lower())
                if matched is None:
                    available = list(mat_id_to_name.values())
                    errors.append(
                        f"{label}: could not resolve material_name {name!r} to an ID. "
                        f"Available materials: {available}"
                    )
                    continue
                s["material_id"] = matched
            else:
                s["material_id"] = int(s["material_id"])
            # Fill in name if missing
            if "material_name" not in s:
                s["material_name"] = mat_id_to_name.get(s["material_id"], f"Material_{s['material_id']}")

        if ptype in ("inlet_q", "exit_h") and "bc_id" not in s:
            errors.append(f"{label}: 'bc_id' is required for type={ptype!r}.")
            continue

        # Validate bounds
        lo = s.get("min")
        hi = s.get("max")
        if lo is None or hi is None:
            errors.append(f"{label}: 'min' and 'max' are required.")
            continue
        lo, hi = float(lo), float(hi)
        if lo >= hi:
            errors.append(f"{label}: 'min' ({lo}) must be less than 'max' ({hi}).")
            continue
        s["min"] = lo
        s["max"] = hi

        # Fill defaults
        midpoint = (lo + hi) / 2.0
        s.setdefault("initial", round(midpoint, 6))
        # Monte Carlo distribution fields are only added when the user provides
        # at least one of them — they are not needed for calibration.
        if any(k in s for k in ("distribution", "mean", "std")):
            s.setdefault("distribution", "truncated_normal")
            s.setdefault("mean", round(midpoint, 6))
            s.setdefault("std", round((hi - lo) / 6.0, 6))

        result_specs.append(s)

    if errors:
        return _err("param_specs validation failed:\n" + "\n".join(errors))

    summary = []
    for s in result_specs:
        if s["type"] == "manning_n":
            summary.append(
                f"manning_n: {s['material_name']} (id={s['material_id']}) "
                f"[{s['min']}, {s['max']}] initial={s['initial']}"
            )
        else:
            summary.append(
                f"{s['type']}: bc_id={s['bc_id']} [{s['min']}, {s['max']}] initial={s['initial']}"
            )

    return _ok(
        result_specs,
        f"Built {len(result_specs)} param spec(s):\n" + "\n".join(f"  {l}" for l in summary),
    )


# ── config tool ───────────────────────────────────────────────────────────────

def load_config(config_file: str = "hmt_config.json") -> dict:
    """Load solver paths and project defaults from a JSON config file.

    This is the recommended way to tell the AI agent about machine-specific
    installation paths (SRH-2D or HEC-RAS executables) and project-level
    defaults (base_case directory, work directory) without repeating them in
    every tool call.

    Create the config file once per machine/project and point the agent to it.
    All values are stored in the session and inherited automatically by
    run_preprocessing(), run_simulation(), run_calibration(), and
    run_monte_carlo().

    Config file format (JSON)
    -------------------------
    .. code-block:: json

        {
          "srh_pre_path":   "C:/Program Files/SMS 13.4 64-bit/python/Lib/site-packages/srh2d_exe/SRH_Pre_Console.exe",
          "srh_path":       "C:/Program Files/SMS 13.4 64-bit/python/Lib/site-packages/srh2d_exe/SRH-2D_Console.exe",
          "extra_dll_path": "C:/Program Files/SMS 13.4 64-bit/python/Lib/site-packages/srh2d_exe",
          "srh_version":    "3.3",
          "hecras_version": "6.6",
          "base_case_dir":  "./base_case",
          "work_dir":       "./runs"
        }

    All keys are optional. Unrecognised keys are ignored.

    Parameters
    ----------
    config_file : str
        Path to the JSON config file. Defaults to ``hmt_config.json`` in the
        current working directory.

    Returns
    -------
    dict with ``"data"`` containing the key/value pairs that were loaded and
    stored in the session.
    """
    config_path = os.path.abspath(config_file)
    if not os.path.isfile(config_path):
        return _err(
            f"Config file not found: {config_path}\n"
            "Create a JSON file with keys such as 'srh_pre_path', 'srh_path', "
            "'extra_dll_path', 'base_case_dir', 'work_dir'. All keys are optional."
        )

    try:
        with open(config_path) as fh:
            cfg = json.load(fh)
    except Exception as exc:
        return _err(f"Failed to parse config file: {exc}")

    session = get_session()
    loaded = {}

    _PATH_KEYS = {
        "srh_pre_path":   "_srh_pre_path",
        "srh_path":       "_srh_path",
        "extra_dll_path": "_extra_dll_path",
        "srh_version":    "_srh_version",
        "hecras_version": "_hecras_version",
        "base_case_dir":  "_base_case_dir",
        "work_dir":       "_work_dir",
    }

    for cfg_key, session_attr in _PATH_KEYS.items():
        if cfg_key in cfg:
            val = cfg[cfg_key]
            setattr(session, session_attr, val)
            loaded[cfg_key] = val

    if not loaded:
        return _err("Config file contained no recognised keys. "
                    f"Valid keys: {list(_PATH_KEYS)}")

    return _ok(
        loaded,
        f"Loaded {len(loaded)} setting(s) from {config_path}: "
        + ", ".join(f"{k}={v!r}" for k, v in loaded.items()),
    )
