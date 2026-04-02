# -*- coding: utf-8 -*-
"""
AI Tools Example — SRH-2D — 02: Modify Parameters and Run

Demonstrates how an AI agent can:
  1. Load solver paths from a config file.
  2. Open a project and check current Manning's n values.
  3. Modify Manning's n for selected material zones.
  4. Modify inlet flow rate.
  5. Save the modified input files.
  6. Run preprocessing and simulation.
  7. Check simulation status.

Edit ../hmt_config.json once to set your SRH-2D installation paths.
Requires SRH-2D to be installed on the machine.
"""

from pyHMT2D.AI_Tools.tools import (
    load_config,
    get_project_info,
    get_materials,
    get_boundary_conditions,
    set_manning_n,
    set_inlet_flow,
    save_modified_inputs,
    run_preprocessing,
    run_simulation,
    get_simulation_status,
)
import json

# ── Configuration ─────────────────────────────────────────────────────────────
CONFIG_FILE  = r"..\hmt_config.json"
CONTROL_FILE = r".\Muncie.srhhydro"
# ─────────────────────────────────────────────────────────────────────────────


def pretty(r): print(json.dumps(r, indent=2, default=str))


# Load config — solver paths (srh_path, srh_pre_path, extra_dll_path) are
# stored in the session and inherited by all subsequent tool calls.
print("=== Load config ===")
pretty(load_config(CONFIG_FILE))

# Open project (sets session model_type = "SRH-2D")
pretty(get_project_info(CONTROL_FILE))

# Check current Manning's n values and boundary conditions
print("\nCurrent materials:")
pretty(get_materials())

print("\nCurrent boundary conditions:")
pretty(get_boundary_conditions())

# Modify Manning's n for materials 2 and 4
print("\nSetting Manning's n ...")
pretty(set_manning_n(
    material_ids=[2, 4],
    values=[0.035, 0.055],
    material_names=["channel", "left_2"],
))

# Modify inlet flow (boundary ID 2)
print("\nSetting inlet flow ...")
pretty(set_inlet_flow(bc_ids=[2], values=[1250.0]))

# Save modified input files
print("\nSaving inputs ...")
pretty(save_modified_inputs())

# Run pre-processor (inherits srh_pre_path and extra_dll_path from session)
print("\nRunning pre-processing ...")
pretty(run_preprocessing(version="3.7"))

# Run solver (inherits srh_path, srh_pre_path, extra_dll_path from session)
print("\nRunning simulation ...")
pretty(run_simulation(timeout_minutes=60))

# Check status
print("\nSimulation status:")
pretty(get_simulation_status())
