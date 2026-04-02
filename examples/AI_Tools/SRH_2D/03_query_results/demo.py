# -*- coding: utf-8 -*-
"""
AI Tools Example — SRH-2D — 03: Query Simulation Results

Demonstrates how an AI agent can interrogate results after a simulation:
  1. Load config and open the project.
  2. Load results from an XMDF file.
  3. Get domain-wide statistics (min/max/mean WSE, depth).
  4. Get WSE at a specific (x, y) location.
  5. Compute the flood extent above a depth threshold.
  6. Extract a cross-section profile.
  7. Export results to VTK for ParaView visualisation.

Edit ../hmt_config.json once to set your installation paths.
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
)
import json

# ── Configuration ─────────────────────────────────────────────────────────────
CONFIG_FILE  = r"..\hmt_config.json"
CONTROL_FILE = r".\Muncie.srhhydro"
RESULT_FILE  = r".\Muncie_XMDFC.h5"
VTK_OUTPUT   = r".\vtk_output"
# ─────────────────────────────────────────────────────────────────────────────


def pretty(r): print(json.dumps(r, indent=2, default=str))


# Load config and open project
print("=== Load config ===")
pretty(load_config(CONFIG_FILE))

pretty(get_project_info(CONTROL_FILE))
pretty(read_results(RESULT_FILE))

# Print the list of available result variables
print("\n=== Available result variables ===")
pretty(get_result_variables(RESULT_FILE))

# Domain-wide statistics — use variable names from get_result_variables() above
print("\n=== WSE Statistics ===")
pretty(get_result_statistics("Water_Elev_ft"))

print("\n=== Depth Statistics ===")
pretty(get_result_statistics("Water_Depth_ft"))

# WSE at a specific point (e.g., a bridge pier location)
print("\n=== WSE at bridge pier ===")
pretty(get_value_at_point(x=411100.0, y=1803450.0, variable="Water_Elev_ft"))

# Flood extent — area with depth >= 0.1 ft
print("\n=== Flood extent (depth >= 0.1 ft) ===")
pretty(get_flood_extent(depth_threshold=0.1))

# Cross-section profile across the main channel
print("\n=== Cross-section profile ===")
result = get_cross_section_profile(
    x1=410500.0, y1=1803000.0,
    x2=411500.0, y2=1803000.0,
    variable="Water_Elev_ft",
    n_points=30,
)
if result["status"] == "ok":
    pts = result["data"]
    print(f"  First point: {pts[0]}")
    print(f"  Last  point: {pts[-1]}")
    print(f"  Total points: {len(pts)}")
else:
    print(f"  ERROR: {result['message']}")

# Export to VTK
print("\n=== Export to VTK ===")
pretty(export_to_vtk(timestep=-1, output_dir=VTK_OUTPUT))
