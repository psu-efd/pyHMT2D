# -*- coding: utf-8 -*-
"""
AI Tools Example — SRH-2D — 04: Manning's n Calibration

Demonstrates AI-assisted calibration of Manning's n against high water marks.

Two usage patterns are shown:
  A. Single evaluation at initial parameters with evaluate_parameters().
  B. Full automated calibration with run_calibration() (GP-minimize loop).

Directory layout expected:
  04_calibration/
  ├── base_case/          <- original, unmodified model files
  │   ├── Muncie.srhhydro
  │   ├── Muncie.srhgeom
  │   └── ...
  ├── HWMs.dat            <- observation file (x, y, wse per row)
  └── demo.py

Edit ../hmt_config.json once to set your SRH-2D installation paths.
"""

from pyHMT2D.AI_Tools.tools import (
    load_config,
    get_project_info,
    build_param_specs,
    check_observation_format,
    load_observations,
    evaluate_parameters,
    run_calibration,
)
import json

# ── Configuration ─────────────────────────────────────────────────────────────
CONFIG_FILE   = r"..\hmt_config.json"
BASE_CASE_DIR = r".\base_case"
CONTROL_FILE  = BASE_CASE_DIR + r"\Muncie.srhhydro"
WORK_DIR      = r"."
OBS_CSV       = r".\HWMs.dat"
OUTPUT_CSV    = r".\calibration_history.csv"
# ─────────────────────────────────────────────────────────────────────────────


def pretty(r): print(json.dumps(r, indent=2, default=str))


# Load config — stores srh_path, srh_pre_path, extra_dll_path in session
print("=== Load config ===")
pretty(load_config(CONFIG_FILE))

# Open project (auto-detects SRH-2D)
pretty(get_project_info(CONTROL_FILE))

# Validate observation file format before proceeding
print("\n=== Check observation file format ===")
fmt = check_observation_format(OBS_CSV)
pretty(fmt)
if fmt["status"] != "ok":
    raise SystemExit("Observation file format error — fix the file before running calibration.")

# Preview observations
print("\n=== Observation points ===")
pretty(load_observations(OBS_CSV))

# ── Build param specs ─────────────────────────────────────────────────────────
# build_param_specs validates material names against the open project and fills
# in defaults — the AI agent calls this instead of constructing dicts manually.
print("\n=== Build parameter specifications ===")
param_result = build_param_specs([
    {
        "type": "manning_n",
        "material_name": "channel",
        "initial": 0.03,
        "min": 0.02,
        "max": 0.08,
    },
    {
        "type": "manning_n",
        "material_name": "left_2",
        "initial": 0.03,
        "min": 0.03,
        "max": 0.12,
    },
])
pretty(param_result)
if param_result["status"] != "ok":
    raise SystemExit("Parameter spec error — check material names match the project.")

PARAM_SPECS = param_result["data"]

# ── Pattern A: single evaluation at initial parameters ────────────────────────
print("\n=== Pattern A: single evaluation at initial parameters ===")
test_specs = [dict(s, value=s["initial"]) for s in PARAM_SPECS]
pretty(evaluate_parameters(
    param_specs=test_specs,
    observation_csv=OBS_CSV,
    result_variable="Water_Elev_ft",
    base_case_dir=BASE_CASE_DIR,
    work_dir=WORK_DIR,
))
# result["data"]["rmse"]      → RMSE at the initial guess
# result["data"]["per_point"] → simulated vs measured at each HWM

# ── Pattern B: full automated GP calibration ──────────────────────────────────
# n_iterations sets the computing budget (total number of model runs).
# Omit it to use the auto-computed default: max(20, 10 * n_params + 5) for "gp".
# Set it explicitly to match your available time/compute, e.g. n_iterations=25.
print("\n=== Pattern B: automated GP calibration ===")
pretty(run_calibration(
    param_specs=PARAM_SPECS,
    observation_csv=OBS_CSV,
    n_iterations=11,        # optional — omit to use auto-computed default
    method="gp",
    result_variable="Water_Elev_ft",
    base_case_dir=BASE_CASE_DIR,
    work_dir=WORK_DIR,
    output_csv=OUTPUT_CSV,
))
# result["data"]["best_params"]           → {"manning_n_2": 0.034, "manning_n_4": 0.071}
# result["data"]["best_rmse"]             → 0.087
# result["data"]["n_iterations_planned"]  → 25
# result["data"]["history_csv"]           → path to full iteration log
