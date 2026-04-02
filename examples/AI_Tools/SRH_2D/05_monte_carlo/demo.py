# -*- coding: utf-8 -*-
"""
AI Tools Example — SRH-2D — 05: Monte Carlo Uncertainty Analysis

Demonstrates AI-assisted Monte Carlo simulation:
  1. Load config and open the project.
  2. Build parameter distributions via build_param_specs().
  3. Generate N samples from truncated-normal distributions.
  4. Run all N cases (serial or parallel multiprocessing).
  5. Compute exceedance-probability statistics at probe points and spatially.

Directory layout expected:
  05_monte_carlo/
  ├── base_case/          <- original, unmodified model files
  │   ├── Muncie.srhhydro
  │   └── ...
  └── demo.py

Edit ../hmt_config.json once to set your SRH-2D installation paths.
"""

from pyHMT2D.AI_Tools.tools import (
    load_config,
    get_project_info,
    build_param_specs,
    generate_mc_samples,
    run_monte_carlo,
    get_mc_statistics,
)
import json

# ── Configuration ─────────────────────────────────────────────────────────────
CONFIG_FILE   = r"..\hmt_config.json"
BASE_CASE_DIR = r".\base_case"
CONTROL_FILE  = BASE_CASE_DIR + r"\Muncie.srhhydro"
WORK_DIR      = r".\mc_runs"
N_SAMPLES     = 100
N_PROCESSES   = 1   # set > 1 to use parallel multiprocessing
# ─────────────────────────────────────────────────────────────────────────────


def pretty(r): print(json.dumps(r, indent=2, default=str))


# Load config — stores srh_path, srh_pre_path, extra_dll_path in session
print("=== Load config ===")
pretty(load_config(CONFIG_FILE))

# Open project (auto-detects SRH-2D)
pretty(get_project_info(CONTROL_FILE))

# ── Build parameter distributions ─────────────────────────────────────────────
print("\n=== Build parameter specifications ===")
param_result = build_param_specs([
    {
        "type": "manning_n",
        "material_name": "channel",
        "distribution": "truncated_normal",
        "mean": 0.04,
        "std": 0.005,
        "min": 0.03,
        "max": 0.05,
    },
])
pretty(param_result)
if param_result["status"] != "ok":
    raise SystemExit("Parameter spec error — check material names match the project.")

PARAM_SPECS = param_result["data"]

# Step 1: Generate samples and save to CSV for audit/reproducibility
print("\n=== Step 1: Generate MC samples ===")
sample_result = generate_mc_samples(
    param_specs=PARAM_SPECS,
    n_samples=N_SAMPLES,
    random_seed=42,
    output_csv=WORK_DIR + r"\mc_samples.csv",
)
pretty(sample_result)

# Step 2: Run all cases using the exact samples generated above
print("\n=== Step 2: Run Monte Carlo ensemble ===")
mc_result = run_monte_carlo(
    base_case_dir=BASE_CASE_DIR,
    param_specs=PARAM_SPECS,
    sample_csv=sample_result["data"]["sample_csv"],  # use Step 1 samples directly
    n_processes=N_PROCESSES,
    delete_cases=True,
    output_dir=WORK_DIR,
)
pretty(mc_result)

# Step 3: Compute exceedance probabilities at probe points and spatially
print("\n=== Step 3: Compute MC statistics ===")
pretty(get_mc_statistics(
    results_json=mc_result["data"]["results_json"],
    variable="Water_Depth_ft",
    observation_points=[
        {"name": "bridge_pier", "x": 411100.0, "y": 1803450.0},
    ],
    exceedance_probabilities=[99, 90, 50, 10, 1],
    output_vtk=WORK_DIR + r"\mc_exceedance.vtk",
))
# result["data"]["point_stats"] → depth at each exceedance prob for each point
# result["data"]["spatial_vtk"] → open in ParaView to see spatial CDFs
