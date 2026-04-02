# -*- coding: utf-8 -*-
"""
AI Tools Example — HEC-RAS 2D — 02: Modify Manning's n and Run

Demonstrates how an AI agent modifies HEC-RAS Manning's n values and
runs the model via COM automation.

Note: set_inlet_flow() and set_exit_wse() are SRH-2D only.
      Manning's n for HEC-RAS is modified directly in the .hdf geometry file.

Requires HEC-RAS installed on Windows.
Edit ../hmt_config.json once to set your HEC-RAS version.
"""

from pyHMT2D.AI_Tools.tools import (
    load_config,
    get_project_info,
    get_materials,
    set_manning_n,
    save_modified_inputs,
    run_simulation,
    get_simulation_status,
    exit_model
)
import json


# ── Configuration ─────────────────────────────────────────────────────────────
CONFIG_FILE = r"..\hmt_config.json"
HDF_FILE    = r".\Muncie\Muncie2D.p01.hdf"
# ─────────────────────────────────────────────────────────────────────────────


def pretty(r): print(json.dumps(r, indent=2, default=str))


# Load config — stores hecras_version in session
print("=== Load config ===")
pretty(load_config(CONFIG_FILE))

# Open project — auto-detects HEC-RAS from .hdf extension
pretty(get_project_info(HDF_FILE))

# Inspect current Manning's n zones
print("\nCurrent Manning's n zones:")
pretty(get_materials())

# Modify Manning's n for zone 1 (channel) and zone 3 (floodplain:left_2)
print("\nSetting Manning's n ...")
pretty(set_manning_n(
    material_ids=[1, 3],
    values=[0.031, 0.062],
    material_names=["channel", "left_2"],
))

# For HEC-RAS the HDF is modified in-place; save_modified_inputs confirms.
print("\nConfirming save ...")
pretty(save_modified_inputs())

# Run HEC-RAS via COM interface (inherits hecras_version from session)
print("\nRunning HEC-RAS simulation ...")
pretty(run_simulation(faceless=False, timeout_minutes=60))

print("\nStatus:")
pretty(get_simulation_status())



# Cleanup: exit HEC-RAS if still running
print("\n=== Cleanup: exit_model ===")
pretty(exit_model())
