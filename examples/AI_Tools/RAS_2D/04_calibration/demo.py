# -*- coding: utf-8 -*-
"""
AI Tools Example — HEC-RAS 2D — 04: Manning's n Calibration

Same calibration tools as SRH-2D. The session auto-routes to HEC-RAS code.

Directory layout expected:
  04_calibration/
  ├── base_case/          <- original, unmodified HEC-RAS files
  │   └── Muncie2D.p01.hdf
  ├── HWMs.dat            <- observation file (x, y, wse per row)
  └── demo.py

Edit ../hmt_config.json once to set your HEC-RAS version.
"""

from pyHMT2D.AI_Tools.tools import (
    load_config,
    get_project_info,
    build_param_specs,
    check_observation_format,
    load_observations,
    evaluate_parameters,
    run_calibration,
    exit_model,
)
import json

# ── Configuration ─────────────────────────────────────────────────────────────
CONFIG_FILE   = r"..\hmt_config.json"
BASE_CASE_DIR = r".\base_case"
HDF_FILE      = BASE_CASE_DIR + r"\Muncie2D.p01.hdf"
WORK_DIR      = r"."
OBS_CSV       = r".\HWMs.dat"
OUTPUT_CSV    = r".\calibration_history.csv"
# ─────────────────────────────────────────────────────────────────────────────


def pretty(r): print(json.dumps(r, indent=2, default=str))


# Load config — stores hecras_version in session
print("=== Load config ===")
pretty(load_config(CONFIG_FILE))

# Open project — auto-detects HEC-RAS
pretty(get_project_info(HDF_FILE))

# Validate observation file format
print("\n=== Check observation file format ===")
fmt = check_observation_format(OBS_CSV)
pretty(fmt)
if fmt["status"] != "ok":
    raise SystemExit("Observation file format error — fix the file before running calibration.")

# Preview observations
print("\n=== Observation points ===")
pretty(load_observations(OBS_CSV))

# Build and validate parameter specs against the open project
print("\n=== Build parameter specifications ===")
param_result = build_param_specs([
    {
        "type": "manning_n",
        "material_name": "channel",
        "initial": 0.04,
        "min": 0.02,
        "max": 0.08,
    },
    {
        "type": "manning_n",
        "material_name": "left_2",
        "initial": 0.06,
        "min": 0.03,
        "max": 0.12,
    },
])
pretty(param_result)
if param_result["status"] != "ok":
    raise SystemExit("Parameter spec error — check material names match the project.")

PARAM_SPECS = param_result["data"]

# Single evaluation at initial parameters
print("\n=== Single evaluation at initial guess ===")
test_specs = [dict(s, value=s["initial"]) for s in PARAM_SPECS]
pretty(evaluate_parameters(
    param_specs=test_specs,
    observation_csv=OBS_CSV,
    result_variable="Water_Elev_ft",
    base_case_dir=BASE_CASE_DIR,
    work_dir=WORK_DIR,
))

# Full GP calibration.
# n_iterations sets the computing budget (total number of model runs).
# Omit it to use the auto-computed default: max(20, 10 * n_params + 5) for "gp".
print("\n=== GP calibration ===")
pretty(run_calibration(
    param_specs=PARAM_SPECS,
    observation_csv=OBS_CSV,
    n_iterations=25,        # optional — omit to use auto-computed default
    method="gp",
    result_variable="Water_Elev_ft",
    base_case_dir=BASE_CASE_DIR,
    work_dir=WORK_DIR,
    output_csv=OUTPUT_CSV,
))

# Cleanup: close the HEC-RAS COM instance so HEC-RAS does not stay open
print("\n=== Cleanup: exit_model ===")
pretty(exit_model())
