# -*- coding: utf-8 -*-
"""
AI Tools Example — HEC-RAS 2D — 01: Open and Inspect a Project

Demonstrates discovery tools for a HEC-RAS project:
  1. Load config.
  2. Find HEC-RAS files in a directory.
  3. Open the HDF result file to retrieve plan metadata.
  4. List all Manning's n (land cover) zones.
  5. List boundary conditions.
  6. Inspect result variables and time steps.

Requires at least one completed HEC-RAS run so a .hdf result file exists.
Edit ../hmt_config.json once to set your HEC-RAS version.
"""

from pyHMT2D.AI_Tools.tools import (
    load_config,
    list_model_files,
    get_project_info,
    get_materials,
    get_boundary_conditions,
    get_result_variables,
    exit_model,
)
import json

# ── Configuration ─────────────────────────────────────────────────────────────
CONFIG_FILE  = r"..\hmt_config.json"
PROJECT_DIR  = r".\Muncie"
HDF_FILE     = r".\Muncie\Muncie2D.p01.hdf"
# ─────────────────────────────────────────────────────────────────────────────


def pretty(r): print(json.dumps(r, indent=2, default=str))


# Step 0: Load config
print("=== Step 0: load_config ===")
pretty(load_config(CONFIG_FILE))

# Step 1: Discover model files
print("\n=== Step 1: list_model_files ===")
pretty(list_model_files(PROJECT_DIR))

# Step 2: Open the HDF file — auto-detects HEC-RAS, stores in session
print("\n=== Step 2: get_project_info ===")
pretty(get_project_info(HDF_FILE))

# Step 3: Manning's n zones
print("\n=== Step 3: get_materials ===")
pretty(get_materials())

# Step 4: Boundary conditions
print("\n=== Step 4: get_boundary_conditions ===")
pretty(get_boundary_conditions())

# Step 5: Result variables
print("\n=== Step 5: get_result_variables ===")
pretty(get_result_variables(HDF_FILE))

# Cleanup: close the HEC-RAS COM instance so HEC-RAS does not stay open
print("\n=== Cleanup: exit_model ===")
pretty(exit_model())
