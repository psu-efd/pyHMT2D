# -*- coding: utf-8 -*-
"""
Group 2 — Parameter modification tools.

All tools read the current session to determine the model type so the caller
never has to specify it. Changes are held in memory until save_modified_inputs()
is called.
"""

from __future__ import annotations

import os
from typing import List, Optional

from pyHMT2D.AI_Tools.state import get_session


def _ok(data, message: str) -> dict:
    return {"status": "ok", "data": data, "message": message}


def _err(message: str) -> dict:
    return {"status": "error", "data": None, "message": message}


# ── tool 6 ────────────────────────────────────────────────────────────────────

def set_manning_n(material_ids: List[int], values: List[float],
                  material_names: Optional[List[str]] = None) -> dict:
    """Set Manning's roughness coefficient for one or more material zones.

    Manning's n controls flow resistance. Typical values: 0.02–0.05 (channel),
    0.04–0.15 (floodplain/vegetation). Changes are held in memory; call
    save_modified_inputs() to write them to disk.

    Use get_materials() first to discover valid material IDs and current values.

    Parameters
    ----------
    material_ids : list of int
        Material zone IDs to update (from get_materials()).
    values : list of float
        New Manning's n values, one per ID. Must be positive (typically 0.01–0.15).
    material_names : list of str, optional
        Zone names for confirmation logging. If omitted, generic names are used.

    Returns
    -------
    dict with "data" containing a list of
        {"id": int, "name": str, "old_n": float, "new_n": float}
    """
    session = get_session()
    try:
        session.require_open()
    except RuntimeError as e:
        return _err(str(e))

    if len(material_ids) != len(values):
        return _err(
            f"material_ids ({len(material_ids)}) and values ({len(values)}) "
            "must have the same length."
        )

    for v in values:
        if not (0.001 < float(v) < 1.0):
            return _err(
                f"Manning's n value {v} is outside the physically plausible range "
                "(0.001, 1.0). Please check your inputs."
            )

    if material_names is None:
        material_names = [f"Material_{mid}" for mid in material_ids]

    try:
        # Capture old values for the return payload
        if session.model_type == "SRH-2D":
            n_dict = session.data.get_ManningN_dict()
            old_values = [float(n_dict.get(mid, float("nan"))) for mid in material_ids]
        else:  # HEC-RAS: ManningNZones = {id: [name_bytes, n_value]}
            zones = session.data.ManningNZones or {}
            old_values = [float(zones[mid][1]) if mid in zones else float("nan") for mid in material_ids]

        # Apply update
        ids_int = [int(i) for i in material_ids]
        vals_float = [float(v) for v in values]

        if session.model_type == "SRH-2D":
            session.data.modify_ManningsNs(ids_int, vals_float, material_names)
        else:  # HEC-RAS
            session.data.modify_ManningsN(ids_int, vals_float, material_names)

        changes = [
            {
                "id": mid,
                "name": name,
                "old_n": old,
                "new_n": new,
            }
            for mid, name, old, new in zip(material_ids, material_names, old_values, vals_float)
        ]

        return _ok(
            changes,
            f"Updated Manning's n for {len(material_ids)} material zone(s). "
            "Call save_modified_inputs() to persist to disk.",
        )

    except Exception as exc:
        return _err(f"Failed to set Manning's n: {exc}")


# ── tool 7 ────────────────────────────────────────────────────────────────────

def set_inlet_flow(bc_ids: List[int], values: List[float]) -> dict:
    """Set inlet flow rate (discharge) for one or more INLET-Q boundaries.

    Applies to SRH-2D only (HEC-RAS inlet flows are set in plan files outside
    pyHMT2D's scope). Use get_boundary_conditions() to find valid BC IDs and
    their current values.

    Parameters
    ----------
    bc_ids : list of int
        Boundary condition IDs with type INLET-Q (from get_boundary_conditions()).
    values : list of float
        New flow rates in m³/s (or ft³/s if the project uses English units).
        Must be non-negative.

    Returns
    -------
    dict with "data" containing a list of
        {"id": int, "old_value": float, "new_value": float, "units": "m3/s"}
    """
    session = get_session()
    try:
        session.require_open()
    except RuntimeError as e:
        return _err(str(e))

    if session.model_type != "SRH-2D":
        return _err(
            "set_inlet_flow() is only supported for SRH-2D. "
            "For HEC-RAS, modify the flow file (.f*) through HEC-RAS directly."
        )

    if len(bc_ids) != len(values):
        return _err(
            f"bc_ids ({len(bc_ids)}) and values ({len(values)}) must have the same length."
        )

    for v in values:
        if float(v) < 0:
            return _err(f"Flow rate {v} must be non-negative.")

    try:
        ids_int   = [int(i)   for i in bc_ids]
        vals_float = [float(v) for v in values]

        session.data.modify_InletQ(ids_int, vals_float)

        changes = [
            {"id": bid, "new_value": val, "units": "m3/s"}
            for bid, val in zip(ids_int, vals_float)
        ]
        return _ok(
            changes,
            f"Updated inlet flow for {len(bc_ids)} boundary(ies). "
            "Call save_modified_inputs() to persist to disk.",
        )

    except Exception as exc:
        return _err(f"Failed to set inlet flow: {exc}")


# ── tool 8 ────────────────────────────────────────────────────────────────────

def set_exit_wse(bc_ids: List[int], values: List[float]) -> dict:
    """Set exit water surface elevation for one or more EXIT-H boundaries.

    Applies to SRH-2D only. Use get_boundary_conditions() to find valid
    EXIT-H boundary IDs and their current values.

    Parameters
    ----------
    bc_ids : list of int
        Boundary condition IDs with type EXIT-H (from get_boundary_conditions()).
    values : list of float
        New water surface elevations in metres (or feet for English-unit projects).

    Returns
    -------
    dict with "data" containing a list of
        {"id": int, "old_value": float, "new_value": float, "units": "m"}
    """
    session = get_session()
    try:
        session.require_open()
    except RuntimeError as e:
        return _err(str(e))

    if session.model_type != "SRH-2D":
        return _err(
            "set_exit_wse() is only supported for SRH-2D. "
            "For HEC-RAS, modify the boundary condition in the plan file."
        )

    if len(bc_ids) != len(values):
        return _err(
            f"bc_ids ({len(bc_ids)}) and values ({len(values)}) must have the same length."
        )

    try:
        ctrl = (
            session.data.srhhydro_obj
            if session.data.control_type == "SRHHydro"
            else session.data.srhsif_obj
        )

        if session.data.control_type == "SRHHydro":
            exh_dict = ctrl.get_ExitH()
        else:
            exh_dict = ctrl.srhsif_content.get("EWSParamsC", {})

        old_values = []
        for bid in bc_ids:
            entry = exh_dict.get(bid)
            if entry is None:
                old_values.append(float("nan"))
            elif isinstance(entry, dict):
                old_values.append(float(entry.get("wse", float("nan"))))
            else:
                old_values.append(float(entry))

        ids_int = [int(i) for i in bc_ids]
        vals_float = [float(v) for v in values]
        session.data.modify_ExitH(ids_int, vals_float)

        changes = [
            {"id": bid, "old_value": old, "new_value": new, "units": "m"}
            for bid, old, new in zip(bc_ids, old_values, vals_float)
        ]

        return _ok(
            changes,
            f"Updated exit WSE for {len(bc_ids)} boundary(ies). "
            "Call save_modified_inputs() to persist to disk.",
        )

    except Exception as exc:
        return _err(f"Failed to set exit WSE: {exc}")


# ── tool 9 ────────────────────────────────────────────────────────────────────

def save_modified_inputs(output_path: Optional[str] = None) -> dict:
    """Write any in-memory parameter changes back to the model input files.

    For SRH-2D: rewrites the .srhhydro or _SIF.dat control file.
    For HEC-RAS: writes the modified Manning's n back to the .hdf geometry.

    Parameters
    ----------
    output_path : str, optional
        If provided, saves to this path instead of overwriting the original.
        For SRH-2D this is the new .srhhydro / _SIF.dat filename.
        For HEC-RAS the HDF is always saved in-place (output_path is ignored).

    Returns
    -------
    dict with "data" containing {"saved_to": str}
    """
    session = get_session()
    try:
        session.require_open()
    except RuntimeError as e:
        return _err(str(e))

    try:
        if session.model_type == "SRH-2D":
            ctrl = (
                session.data.srhhydro_obj
                if session.data.control_type == "SRHHydro"
                else session.data.srhsif_obj
            )
            ctrl.save_as(output_path)
            saved_to = output_path or session.project_file
        else:  # HEC-RAS — changes are already written into the HDF in-place
            # by modify_ManningsN(); just confirm the file path
            saved_to = session.project_file

        return _ok(
            {"saved_to": saved_to},
            f"Input file saved to: {saved_to}",
        )

    except Exception as exc:
        return _err(f"Failed to save inputs: {exc}")
