# -*- coding: utf-8 -*-
"""Tests for AI_Tools monte_carlo_tools (no solver required)."""

import os
import json
import pytest

from pyHMT2D.AI_Tools.state import reset_session
from pyHMT2D.AI_Tools.tools.monte_carlo_tools import (
    generate_mc_samples,
    run_monte_carlo,
    get_mc_statistics,
)


@pytest.fixture(autouse=True)
def clean_session():
    reset_session()
    yield
    reset_session()


# ── generate_mc_samples ───────────────────────────────────────────────────────

def test_generate_mc_samples_basic(tmp_path):
    specs = [
        {
            "type": "manning_n",
            "material_id": 2,
            "distribution": "truncated_normal",
            "mean": 0.04,
            "std": 0.005,
            "min": 0.03,
            "max": 0.05,
        }
    ]
    result = generate_mc_samples(
        specs, n_samples=10, random_seed=0,
        output_csv=str(tmp_path / "samples.csv"),
    )
    assert result["status"] == "ok"
    assert len(result["data"]["samples"]) == 10
    assert len(result["data"]["samples"][0]) == 1  # 1 parameter


def test_generate_mc_samples_uniform(tmp_path):
    specs = [
        {"type": "inlet_q", "bc_id": 1, "distribution": "uniform",
         "min": 50.0, "max": 200.0}
    ]
    result = generate_mc_samples(
        specs, n_samples=5, output_csv=str(tmp_path / "samples.csv")
    )
    assert result["status"] == "ok"
    for row in result["data"]["samples"]:
        assert 50.0 <= row[0] <= 200.0


def test_generate_mc_samples_bounds_respected(tmp_path):
    specs = [
        {"type": "manning_n", "material_id": 1, "distribution": "truncated_normal",
         "mean": 0.04, "std": 0.01, "min": 0.02, "max": 0.06}
    ]
    result = generate_mc_samples(
        specs, n_samples=50, random_seed=99,
        output_csv=str(tmp_path / "s.csv"),
    )
    assert result["status"] == "ok"
    for row in result["data"]["samples"]:
        assert 0.02 <= row[0] <= 0.06


def test_generate_mc_samples_reproducible(tmp_path):
    specs = [{"type": "manning_n", "material_id": 1, "distribution": "uniform",
              "min": 0.02, "max": 0.08}]
    r1 = generate_mc_samples(specs, n_samples=5, random_seed=42,
                              output_csv=str(tmp_path / "a.csv"))
    r2 = generate_mc_samples(specs, n_samples=5, random_seed=42,
                              output_csv=str(tmp_path / "b.csv"))
    assert r1["data"]["samples"] == r2["data"]["samples"]


def test_generate_mc_samples_zero():
    result = generate_mc_samples([], n_samples=0)
    assert result["status"] == "error"


def test_generate_mc_samples_csv_written(tmp_path):
    specs = [{"type": "manning_n", "material_id": 1, "distribution": "uniform",
              "min": 0.02, "max": 0.08}]
    csv_path = str(tmp_path / "out.csv")
    generate_mc_samples(specs, n_samples=5, output_csv=csv_path)
    assert os.path.isfile(csv_path)


# ── run_monte_carlo ───────────────────────────────────────────────────────────

def test_run_monte_carlo_no_project():
    result = run_monte_carlo(
        "/some/base_case",
        [{"type": "manning_n", "material_id": 1,
          "distribution": "uniform", "min": 0.02, "max": 0.08}],
        n_samples=2,
    )
    assert result["status"] == "error"
    assert "no project" in result["message"].lower()


def test_run_monte_carlo_missing_base_case_dir():
    from pyHMT2D.AI_Tools.state import get_session
    s = get_session()
    s.model_type = "SRH-2D"
    s.data = object()
    s.project_file = "C:/fake/Muncie.srhhydro"

    result = run_monte_carlo(
        "/nonexistent/base_case",
        [{"type": "manning_n", "material_id": 1,
          "distribution": "uniform", "min": 0.02, "max": 0.08}],
        n_samples=2,
    )
    assert result["status"] == "error"
    assert "not found" in result["message"].lower()


# ── get_mc_statistics ─────────────────────────────────────────────────────────

def test_get_mc_statistics_missing_json():
    result = get_mc_statistics("/no/such/results.json")
    assert result["status"] == "error"
    assert "not found" in result["message"].lower()


def test_get_mc_statistics_no_vtk_files(tmp_path):
    # Write a results JSON with empty vtk_files list
    rjson = tmp_path / "results.json"
    rjson.write_text(json.dumps({"vtk_files": [], "n_successful": 0}))
    result = get_mc_statistics(str(rjson))
    assert result["status"] == "error"
    assert "no vtk" in result["message"].lower()
