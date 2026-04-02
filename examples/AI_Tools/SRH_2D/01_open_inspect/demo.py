# -*- coding: utf-8 -*-
"""
AI Tools Example — SRH-2D — 01: Open and Inspect a Project

Demonstrates how an AI agent uses the discovery tools to:
  1. Load machine-specific paths from a config file.
  2. Find all model files in a directory.
  3. Open a project and retrieve its metadata.
  4. List all material zones with their Manning's n values.
  5. List all boundary conditions.
  6. Inspect a result file's available variables and time steps.

No solver installation is required for steps 2–6.
Edit ../hmt_config.json once to set your installation paths.
"""

from pyHMT2D.AI_Tools.tools import (
    load_config,
    list_model_files,
    get_project_info,
    get_materials,
    get_boundary_conditions,
    get_result_variables,
)
import json

# ── Configuration ─────────────────────────────────────────────────────────────
# Machine-specific paths (solver executables) are stored in hmt_config.json.
# Only project-specific paths are set here.
CONFIG_FILE  = r"..\hmt_config.json"
PROJECT_DIR  = r"."
CONTROL_FILE = r"Muncie.srhhydro"
RESULT_FILE  = r"Muncie_XMDFC.h5"
# ─────────────────────────────────────────────────────────────────────────────


def pretty(result: dict) -> None:
    print(json.dumps(result, indent=2, default=str))


# Step 0: Load config — stores solver paths in session for later tools
print("=== Step 0: load_config ===")
pretty(load_config(CONFIG_FILE))

# Step 1: Find all model files in the project directory
print("\n=== Step 1: list_model_files ===")
pretty(list_model_files(PROJECT_DIR))

# Step 2: Open the project — auto-detects SRH-2D, stores in session
print("\n=== Step 2: get_project_info ===")
pretty(get_project_info(CONTROL_FILE))

# Step 3: List material zones and their Manning's n values
print("\n=== Step 3: get_materials ===")
pretty(get_materials())

# Step 4: List boundary conditions
print("\n=== Step 4: get_boundary_conditions ===")
pretty(get_boundary_conditions())

# Step 5: Inspect result variables and time steps
print("\n=== Step 5: get_result_variables ===")
pretty(get_result_variables(RESULT_FILE))
