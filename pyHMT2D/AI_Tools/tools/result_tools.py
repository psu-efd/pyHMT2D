# -*- coding: utf-8 -*-
"""
Group 4 — Result query tools.

All spatial queries use the VTK representation of the simulation results so
that both SRH-2D and HEC-RAS share the same code path after export.
"""

from __future__ import annotations

import os
import tempfile
from typing import List, Optional

import numpy as np

from pyHMT2D.AI_Tools.state import get_session


def _ok(data, message: str) -> dict:
    return {"status": "ok", "data": data, "message": message}


def _err(message: str) -> dict:
    return {"status": "error", "data": None, "message": message}


def _get_vtk_cell_array_names(reader) -> list:
    """Return a list of cell data array names present in a VTK unstructured grid."""
    ugrid = reader.GetOutput()
    cd = ugrid.GetCellData()
    return [cd.GetArrayName(i) for i in range(cd.GetNumberOfArrays())]


def _resolve_vtk_varname(name: str, available: list) -> str | None:
    """Resolve a user-supplied variable name to an actual VTK cell-array name.

    Tries in order:
    1. Exact match.
    2. Case-insensitive, spaces/hyphens → underscores normalization.
    3. difflib closest match (ratio >= 0.45) for HEC-RAS HDF names like
       "Water Surface" that differ from VTK names like "Water_Elev_ft".

    Returns the matched VTK name, or None if no reasonable match is found.
    """
    if not available:
        return None

    # 1. Exact
    if name in available:
        return name

    def _norm(s: str) -> str:
        return s.lower().replace(" ", "_").replace("-", "_")

    norm_name = _norm(name)

    # 2. Normalised exact
    norm_map = {_norm(a): a for a in available}
    if norm_name in norm_map:
        return norm_map[norm_name]

    # 3. difflib closest match on normalised names
    import difflib
    candidates = list(norm_map.keys())
    matches = difflib.get_close_matches(norm_name, candidates, n=1, cutoff=0.40)
    if matches:
        return norm_map[matches[0]]

    return None


def _ensure_vtk(session, timestep: int = -1) -> tuple:
    """Export results to a temporary VTK file if not already done.

    Returns (vtk_file_path, error_message_or_None).
    """
    if session.vtk_file and os.path.isfile(session.vtk_file):
        return session.vtk_file, None

    # Need to export
    tmp_dir = tempfile.mkdtemp(prefix="pyhmt2d_ai_")
    try:
        if session.model_type == "SRH-2D":
            session.data.outputXMDFDataToVTK(
                bNodal=False,
                timeStep=timestep,
                lastTimeStep=(timestep == -1),
                dir=tmp_dir,
            )
            # Find the produced VTK file
            import glob
            vtk_files = glob.glob(os.path.join(tmp_dir, "*.vtk"))
            if not vtk_files:
                return None, "SRH-2D VTK export produced no files."
            vtk_file = sorted(vtk_files)[-1]
        else:  # HEC-RAS
            session.data.saveHEC_RAS2D_results_to_VTK(
                timeStep=timestep,
                lastTimeStep=(timestep == -1),
                fileNameBase="ras_result",
                dir=tmp_dir,
            )
            import glob
            vtk_files = glob.glob(os.path.join(tmp_dir, "*.vtk"))
            if not vtk_files:
                return None, "HEC-RAS VTK export produced no files."
            vtk_file = sorted(vtk_files)[-1]

        session.vtk_file = vtk_file
        return vtk_file, None

    except Exception as exc:
        return None, f"VTK export failed: {exc}"


# ── tool 13 ───────────────────────────────────────────────────────────────────

def read_results(result_file: str, timestep: int = -1) -> dict:
    """Load simulation results into the session for subsequent queries.

    For SRH-2D: reads an XMDF/HDF5 result file (.h5 or .xmdf).
    For HEC-RAS: the .hdf plan result file is already loaded when the project
    is opened; call this tool to (re)load solution arrays for a specific
    time step.

    Parameters
    ----------
    result_file : str
        Path to the SRH-2D XMDF file or HEC-RAS HDF file.
    timestep : int
        Time step index to load (-1 = last time step). Default -1.

    Returns
    -------
    dict with "data" containing
        {"variables": [...], "n_time_steps": int, "time_values": [...]}
    """
    if not os.path.isfile(result_file):
        return _err(f"Result file not found: {result_file}")

    session = get_session()
    try:
        session.require_open()
    except RuntimeError as e:
        return _err(str(e))

    try:
        if session.model_type == "SRH-2D":
            session.data.readSRHXMDFFile(result_file, bNodal=False)
            time_arr = session.data.xmdfTimeArray_Cell
            var_names = list(session.data.xmdfAllData_Cell.keys()) if session.data.xmdfAllData_Cell else []
            time_values = time_arr.tolist() if time_arr is not None else []

        else:  # HEC-RAS
            session.data.load2DAreaSolutions()
            time_values_raw = session.data.get2DAreaSolutionTimes()
            time_values = [float(t) for t in time_values_raw] if time_values_raw is not None else []
            # Variable names from HDF
            import h5py
            var_names = []
            with h5py.File(result_file, "r") as hf:
                try:
                    area_grp = (
                        hf["Results"]["Unsteady"]["Output"]["Output Blocks"]
                        ["Base Output"]["Unsteady Time Series"]["2D Flow Areas"]
                    )
                    first_area = list(area_grp.keys())[0]
                    var_names = [
                        k for k in area_grp[first_area].keys()
                        if isinstance(area_grp[first_area][k], h5py.Dataset)
                    ]
                except KeyError:
                    pass

        session.result_file = result_file
        session.result_loaded = True
        session.result_time_steps = time_values
        session.result_variables = var_names
        session.vtk_file = None  # invalidate cached VTK

        info = {
            "result_file": result_file,
            "variables": var_names,
            "n_time_steps": len(time_values),
            "time_values": time_values[:20],
            "time_values_truncated": len(time_values) > 20,
        }
        return _ok(info, f"Loaded {len(var_names)} variable(s) over {len(time_values)} time step(s).")

    except Exception as exc:
        return _err(f"Failed to load results: {exc}")


# ── tool 14 ───────────────────────────────────────────────────────────────────

def get_value_at_point(
    x: float,
    y: float,
    variable: str,
    timestep: int = -1,
) -> dict:
    """Interpolate a result variable at a specific (x, y) coordinate.

    Uses VTK Gaussian-kernel interpolation onto the query point. The result
    file must be loaded first (call read_results()).

    Parameters
    ----------
    x : float
        Easting / x-coordinate in the model's coordinate system.
    y : float
        Northing / y-coordinate in the model's coordinate system.
    variable : str
        Name of the result variable (from get_result_variables() or
        read_results()), e.g. "Water_Elev_m", "Velocity_m_p_s".
    timestep : int
        Time step index (-1 = last step). Default -1.

    Returns
    -------
    dict with "data" containing
        {"x": float, "y": float, "variable": str, "value": float, "units": str}
    """
    session = get_session()
    try:
        session.require_results()
    except RuntimeError as e:
        return _err(str(e))

    try:
        import vtk
        from pyHMT2D.Misc import vtkHandler

        vtk_file, err = _ensure_vtk(session, timestep)
        if err:
            return _err(err)

        handler = vtkHandler()
        reader = handler.readVTK_UnstructuredGrid(vtk_file)

        available = _get_vtk_cell_array_names(reader)
        resolved_var = _resolve_vtk_varname(variable, available)
        if resolved_var is None:
            return _err(
                f"Variable '{variable}' not found in result file. "
                f"Available variables: {available}"
            )
        variable = resolved_var

        # Build a single vtkPoints object
        pts = vtk.vtkPoints()
        pts.InsertNextPoint(float(x), float(y), 0.0)

        probed = handler.probeUnstructuredGridVTKOnPoints(pts, reader, variable)
        # probed returns (points_array, values_array, elev_array)
        if probed is None or len(probed) < 2:
            return _err(f"Probing failed for variable '{variable}' at ({x}, {y}).")

        value = float(probed[1][0]) if hasattr(probed[1], "__len__") else float(probed[1])

        return _ok(
            {"x": x, "y": y, "variable": variable, "value": value, "timestep": timestep},
            f"{variable} = {value:.4f} at ({x:.1f}, {y:.1f})",
        )

    except Exception as exc:
        return _err(f"Failed to get value at point: {exc}")


# ── tool 15 ───────────────────────────────────────────────────────────────────

def get_result_statistics(
    variable: str,
    timestep: int = -1,
    depth_threshold: float = None,
) -> dict:
    """Return domain-wide statistics (min, max, mean, std) for a result variable.

    Parameters
    ----------
    variable : str
        Result variable name (e.g. "Water_Depth_m", "Water_Elev_m").
    timestep : int
        Time step index (-1 = last step). Default -1.
    depth_threshold : float, optional
        If provided, only include cells where water depth >= this value (wet cells).
        The depth variable is auto-detected from the loaded results.

    Returns
    -------
    dict with "data" containing
        {"variable": str, "min": float, "max": float, "mean": float,
         "std": float, "n_cells": int}
    """
    session = get_session()
    try:
        session.require_results()
    except RuntimeError as e:
        return _err(str(e))

    try:
        idx = timestep if timestep >= 0 else -1

        def _load_array(varname):
            """Load a variable array from in-memory data or VTK fallback."""
            if session.model_type == "SRH-2D" and session.data.xmdfAllData_Cell:
                raw = session.data.xmdfAllData_Cell.get(varname)
                if raw is not None:
                    return np.array(raw[idx]).flatten(), None
            from pyHMT2D.Misc import vtkHandler
            vtk_file, err = _ensure_vtk(session, timestep)
            if err:
                return None, err
            handler = vtkHandler()
            reader = handler.readVTK_UnstructuredGrid(vtk_file)
            available = _get_vtk_cell_array_names(reader)
            resolved = _resolve_vtk_varname(varname, available)
            if resolved is None:
                return None, (
                    f"Variable '{varname}' not found in result file. "
                    f"Available variables: {available}"
                )
            arr_vtk = handler.get_uGRid_cell_field_with_name(reader, resolved)
            if arr_vtk is None:
                return None, f"Failed to read array for variable '{resolved}'."
            return np.array(arr_vtk).flatten(), None

        arr, err = _load_array(variable)
        if err:
            return _err(err)

        # Build wet-cell mask if depth_threshold is given
        wet_mask = None
        if depth_threshold is not None:
            # Auto-detect depth variable name
            depth_var = None
            if session.model_type == "SRH-2D" and session.data.xmdfAllData_Cell:
                candidates = list(session.data.xmdfAllData_Cell.keys())
            else:
                candidates = []
            for candidate in candidates:
                if "depth" in candidate.lower():
                    depth_var = candidate
                    break
            if depth_var is None:
                return _err(
                    "Could not auto-detect depth variable for wet-cell masking. "
                    "Available variables: " + str(candidates)
                )
            depth_arr, err = _load_array(depth_var)
            if err:
                return _err(f"Failed to load depth variable '{depth_var}': {err}")
            depth_arr = depth_arr[np.isfinite(depth_arr)]
            wet_mask = depth_arr >= depth_threshold

        # Remove NaN / fill values
        arr = arr[np.isfinite(arr)]
        if len(arr) == 0:
            return _err(f"No valid data for variable '{variable}'.")

        # Apply wet-cell mask
        if wet_mask is not None:
            if len(wet_mask) != len(arr):
                return _err(
                    f"Depth mask length ({len(wet_mask)}) does not match "
                    f"variable array length ({len(arr)})."
                )
            arr = arr[wet_mask]
            if len(arr) == 0:
                return _err(
                    f"No wet cells found with depth >= {depth_threshold}."
                )

        stats = {
            "variable": variable,
            "min": float(np.min(arr)),
            "max": float(np.max(arr)),
            "mean": float(np.mean(arr)),
            "std": float(np.std(arr)),
            "n_cells": int(len(arr)),
            "timestep": timestep,
            "depth_threshold": depth_threshold,
        }
        suffix = (
            f" (wet cells, depth >= {depth_threshold})" if depth_threshold is not None else ""
        )
        return _ok(
            stats,
            f"{variable}: min={stats['min']:.4f}, max={stats['max']:.4f}, "
            f"mean={stats['mean']:.4f} over {stats['n_cells']} cells{suffix}.",
        )

    except Exception as exc:
        return _err(f"Failed to compute statistics: {exc}")


# ── tool 16 ───────────────────────────────────────────────────────────────────

def get_flood_extent(
    depth_threshold: float = 0.01,
    timestep: int = -1,
) -> dict:
    """Calculate the flooded area based on a minimum water depth threshold.

    Sums cell areas where water depth exceeds the threshold. Requires results
    to be loaded (call read_results()).

    Parameters
    ----------
    depth_threshold : float
        Minimum water depth (m) to be considered flooded. Default 0.01 m.
    timestep : int
        Time step index (-1 = last step). Default -1.

    Returns
    -------
    dict with "data" containing
        {"flooded_area_m2": float, "flooded_area_km2": float,
         "total_domain_area_m2": float, "flooded_fraction": float,
         "depth_threshold_m": float, "n_flooded_cells": int}
    """
    session = get_session()
    try:
        session.require_results()
    except RuntimeError as e:
        return _err(str(e))

    try:
        from pyHMT2D.Misc import vtkHandler
        import vtk
        from vtk.util import numpy_support as VN

        vtk_file, err = _ensure_vtk(session, timestep)
        if err:
            return _err(err)

        handler = vtkHandler()
        reader = handler.readVTK_UnstructuredGrid(vtk_file)
        ugrid = reader.GetOutput()

        # Identify depth variable name
        depth_var = None
        for candidate in ["Water_Depth_m", "Depth_m", "Water_Depth_ft", "water_depth", "Depth"]:
            if ugrid.GetCellData().GetArray(candidate) is not None:
                depth_var = candidate
                break
        if depth_var is None:
            # Try to derive from WSE - bed
            depth_var = handler.get_uGRid_cell_field_with_name(reader, "Water_Depth_m")

        depth_arr = handler.get_uGRid_cell_field_with_name(reader, depth_var or "Water_Depth_m")
        if depth_arr is None:
            return _err(
                "Cannot find a water depth variable. "
                "Ensure results include 'Water_Depth_m' or similar."
            )
        depths = np.array(depth_arr).flatten()

        # Compute cell areas using VTK
        n_cells = ugrid.GetNumberOfCells()
        cell_areas = np.zeros(n_cells)
        for ci in range(n_cells):
            cell = ugrid.GetCell(ci)
            pts = cell.GetPoints()
            if pts.GetNumberOfPoints() >= 3:
                p0 = np.array(pts.GetPoint(0)[:2])
                area = 0.0
                for pi in range(1, pts.GetNumberOfPoints() - 1):
                    p1 = np.array(pts.GetPoint(pi)[:2])
                    p2 = np.array(pts.GetPoint(pi + 1)[:2])
                    area += 0.5 * abs(np.cross(p1 - p0, p2 - p0))
                cell_areas[ci] = area

        flooded_mask = depths >= depth_threshold
        flooded_area = float(np.sum(cell_areas[flooded_mask]))
        total_area = float(np.sum(cell_areas))
        fraction = flooded_area / total_area if total_area > 0 else 0.0

        result = {
            "flooded_area_m2": flooded_area,
            "flooded_area_km2": flooded_area / 1e6,
            "total_domain_area_m2": total_area,
            "total_domain_area_km2": total_area / 1e6,
            "flooded_fraction": round(fraction, 4),
            "flooded_percent": round(fraction * 100, 2),
            "depth_threshold_m": depth_threshold,
            "n_flooded_cells": int(np.sum(flooded_mask)),
            "n_total_cells": n_cells,
        }

        return _ok(
            result,
            f"Flooded area: {flooded_area/1e6:.3f} km² "
            f"({fraction*100:.1f}% of domain) at depth ≥ {depth_threshold} m.",
        )

    except Exception as exc:
        return _err(f"Failed to compute flood extent: {exc}")


# ── tool 17 ───────────────────────────────────────────────────────────────────

def get_cross_section_profile(
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    variable: str,
    n_points: int = 50,
    timestep: int = -1,
) -> dict:
    """Extract a result variable along a straight-line transect.

    Interpolates the requested variable at evenly spaced points between
    (x1, y1) and (x2, y2).

    Parameters
    ----------
    x1, y1 : float
        Start coordinates of the transect.
    x2, y2 : float
        End coordinates of the transect.
    variable : str
        Result variable name to extract along the line.
    n_points : int
        Number of sample points along the transect. Default 50.
    timestep : int
        Time step index (-1 = last step). Default -1.

    Returns
    -------
    dict with "data" containing a list of
        {"distance_m": float, "x": float, "y": float, "value": float}
    """
    session = get_session()
    try:
        session.require_results()
    except RuntimeError as e:
        return _err(str(e))

    if n_points < 2:
        return _err("n_points must be at least 2.")

    try:
        import vtk
        from pyHMT2D.Misc import vtkHandler

        vtk_file, err = _ensure_vtk(session, timestep)
        if err:
            return _err(err)

        handler = vtkHandler()
        reader = handler.readVTK_UnstructuredGrid(vtk_file)

        available = _get_vtk_cell_array_names(reader)
        resolved_var = _resolve_vtk_varname(variable, available)
        if resolved_var is None:
            return _err(
                f"Variable '{variable}' not found in result file. "
                f"Available variables: {available}"
            )
        variable = resolved_var

        # Build evenly spaced points
        xs = np.linspace(x1, x2, n_points)
        ys = np.linspace(y1, y2, n_points)
        total_dist = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        distances = np.linspace(0, total_dist, n_points)

        pts = vtk.vtkPoints()
        for xi, yi in zip(xs, ys):
            pts.InsertNextPoint(float(xi), float(yi), 0.0)

        probed = handler.probeUnstructuredGridVTKOnPoints(pts, reader, variable)
        if probed is None or len(probed) < 2:
            return _err(f"Probing failed for variable '{variable}'.")

        values = probed[1]

        profile = [
            {
                "distance_m": round(float(d), 2),
                "x": round(float(xi), 2),
                "y": round(float(yi), 2),
                "value": round(float(v), 4),
            }
            for d, xi, yi, v in zip(distances, xs, ys, values)
        ]

        return _ok(
            profile,
            f"Extracted {variable} at {n_points} points over {total_dist:.1f} m transect.",
        )

    except Exception as exc:
        return _err(f"Failed to extract cross-section profile: {exc}")
