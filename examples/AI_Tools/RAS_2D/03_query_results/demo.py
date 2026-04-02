# -*- coding: utf-8 -*-
"""
AI Tools Example — HEC-RAS 2D — 03: Query HEC-RAS Results

The same result-query tools work for HEC-RAS as for SRH-2D — the session
dispatches to the right code path based on the detected model type.

Edit ../hmt_config.json once to set your HEC-RAS version.
"""

from pyHMT2D.AI_Tools.tools import (
    load_config,
    get_project_info,
    read_results,
    get_result_variables,
    get_result_statistics,
    get_value_at_point,
    get_flood_extent,
    get_cross_section_profile,
    export_to_vtk,
    exit_model
)
import json

# ── Configuration ─────────────────────────────────────────────────────────────
CONFIG_FILE = r"..\hmt_config.json"
HDF_FILE    = r".\Muncie\Muncie2D.p01.hdf"
VTK_OUTPUT  = r".\vtk_output"
# ─────────────────────────────────────────────────────────────────────────────


def pretty(r): print(json.dumps(r, indent=2, default=str))


# Load config and open project
print("=== Load config ===")
pretty(load_config(CONFIG_FILE))

pretty(get_project_info(HDF_FILE))
pretty(read_results(HDF_FILE))

# Check available variable names first — use these names in calls below
print("\n=== Available result variables ===")
pretty(get_result_variables(HDF_FILE))

# Statistics — use variable names from get_result_variables() above
print("\n=== WSE Statistics ===")
pretty(get_result_statistics("Water Surface"))

print("\n=== Depth Statistics ===")
pretty(get_result_statistics("Depth"))

# Point query
print("\n=== WSE at probe point ===")
pretty(get_value_at_point(x=411100.0, y=1803450.0, variable="Water Surface"))

# Flood extent
print("\n=== Flood extent (depth >= 0.1 m) ===")
pretty(get_flood_extent(depth_threshold=0.1))

# Cross-section
print("\n=== Cross-section profile ===")
result = get_cross_section_profile(
    x1=410500.0, y1=1803000.0,
    x2=411500.0, y2=1803000.0,
    variable="Water Surface",
    n_points=30,
)
if result["status"] == "ok":
    pts = result["data"]
    print(f"  Profile has {len(pts)} points; "
          f"min={min(p['value'] for p in pts):.3f}, "
          f"max={max(p['value'] for p in pts):.3f}")
else:
    print(f"  ERROR: {result['message']}")

# Export VTK
print("\n=== Export to VTK ===")
pretty(export_to_vtk(timestep=-1, output_dir=VTK_OUTPUT))

# Cleanup: exit HEC-RAS COM instance so HEC-RAS does not stay open
print("\n=== Cleanup: exit_model ===")
pretty(exit_model())
