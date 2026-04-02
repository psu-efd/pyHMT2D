# -*- coding: utf-8 -*-
"""
AI Tools Example — HEC-RAS 2D — 05: Monte Carlo Uncertainty Analysis

Same MC tools as SRH-2D. The session routes to HEC-RAS execution paths.

Directory layout expected:
  05_monte_carlo/
  ├── base_case/          <- original, unmodified HEC-RAS files
  │   └── Muncie2D.p01.hdf
  └── demo.py

Edit ../hmt_config.json once to set your HEC-RAS version.
"""

from pyHMT2D.AI_Tools.tools import (
    load_config,
    get_project_info,
    build_param_specs,
    generate_mc_samples,
    run_monte_carlo,
    get_mc_statistics,
    exit_model,
)
import json
import os

# ── Configuration ─────────────────────────────────────────────────────────────
CONFIG_FILE   = r"..\hmt_config.json"
BASE_CASE_DIR = r".\base_case"
HDF_FILE      = BASE_CASE_DIR + r"\Muncie2D.p01.hdf"
WORK_DIR      = r".\mc_runs"
N_SAMPLES     = 100
N_PROCESSES   = 1   # set > 1 to use parallel multiprocessing
# ─────────────────────────────────────────────────────────────────────────────


def pretty(r): print(json.dumps(r, indent=2, default=str))


# Load config — stores hecras_version in session
print("=== Load config ===")
pretty(load_config(CONFIG_FILE))

# Open project — auto-detects HEC-RAS
pretty(get_project_info(HDF_FILE))

# Create the mc_runs directory if it doesn't exist
if not os.path.exists(WORK_DIR):
    os.makedirs(WORK_DIR)

# Build and validate parameter distributions
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
    sample_csv=sample_result["data"]["sample_csv"],
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
        {"name": "probe_A", "x": 411100.0, "y": 1803450.0},
    ],
    exceedance_probabilities=[99, 90, 50, 10, 1],
    output_vtk=WORK_DIR + r"\mc_exceedance.vtk",
))

# Cleanup: close the HEC-RAS COM instance so HEC-RAS does not stay open
print("\n=== Cleanup: exit_model ===")
pretty(exit_model())
