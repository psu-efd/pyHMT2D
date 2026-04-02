# -*- coding: utf-8 -*-
"""Tests for AI_Tools calibration_tools (no solver required)."""

import os
import pytest

from pyHMT2D.AI_Tools.state import reset_session
from pyHMT2D.AI_Tools.tools.calibration_tools import (
    load_observations,
    evaluate_parameters,
    run_calibration,
)


@pytest.fixture(autouse=True)
def clean_session():
    reset_session()
    yield
    reset_session()


# ── load_observations ─────────────────────────────────────────────────────────

def test_load_observations_missing_file():
    result = load_observations("/no/such/hwm.csv")
    assert result["status"] == "error"
    assert "not found" in result["message"].lower()


def test_load_observations_3col(tmp_path):
    csv = tmp_path / "hwm.csv"
    csv.write_text(
        "# comment\n"
        "413253.87, 1800957.98, 951.46\n"
        "411429.26, 1803300.59, 945.40\n"
    )
    result = load_observations(str(csv))
    assert result["status"] == "ok"
    assert len(result["data"]) == 2
    assert result["data"][0]["x"] == pytest.approx(413253.87)
    assert result["data"][0]["observed_value"] == pytest.approx(951.46)


def test_load_observations_4col(tmp_path):
    csv = tmp_path / "hwm.csv"
    csv.write_text("pt1, 413253.87, 1800957.98, 951.46\n")
    result = load_observations(str(csv))
    assert result["status"] == "ok"
    assert result["data"][0]["name"] == "pt1"


def test_load_observations_empty(tmp_path):
    csv = tmp_path / "hwm.csv"
    csv.write_text("# only comments\n")
    result = load_observations(str(csv))
    assert result["status"] == "error"


# ── evaluate_parameters ───────────────────────────────────────────────────────

def test_evaluate_parameters_no_project(tmp_path):
    csv = tmp_path / "hwm.csv"
    csv.write_text("413000, 1800000, 950.0\n")
    result = evaluate_parameters(
        [{"type": "manning_n", "material_id": 1, "value": 0.04}],
        str(csv),
    )
    assert result["status"] == "error"
    assert "no project" in result["message"].lower()


def test_evaluate_parameters_missing_obs(tmp_path):
    from pyHMT2D.AI_Tools.state import get_session
    s = get_session()
    s.model_type = "SRH-2D"
    s.data = object()
    s.project_file = str(tmp_path / "test.srhhydro")

    result = evaluate_parameters(
        [{"type": "manning_n", "material_id": 1, "value": 0.04}],
        "/no/hwm.csv",
    )
    assert result["status"] == "error"
    assert "not found" in result["message"].lower()


def test_evaluate_parameters_missing_value_key(tmp_path):
    from pyHMT2D.AI_Tools.state import get_session
    csv_file = tmp_path / "hwm.csv"
    csv_file.write_text("413000, 1800000, 950.0\n")
    s = get_session()
    s.model_type = "SRH-2D"
    s.data = object()
    s.project_file = str(tmp_path / "test.srhhydro")

    result = evaluate_parameters(
        [{"type": "manning_n", "material_id": 1}],  # missing "value"
        str(csv_file),
    )
    assert result["status"] == "error"
    assert "value" in result["message"]


# ── run_calibration ───────────────────────────────────────────────────────────

def test_run_calibration_no_project(tmp_path):
    csv = tmp_path / "hwm.csv"
    csv.write_text("413000, 1800000, 950.0\n")
    result = run_calibration(
        [{"type": "manning_n", "material_id": 1, "initial": 0.04,
          "min": 0.02, "max": 0.08}],
        str(csv),
    )
    assert result["status"] == "error"


def test_run_calibration_bad_method(tmp_path):
    from pyHMT2D.AI_Tools.state import get_session
    csv_file = tmp_path / "hwm.csv"
    csv_file.write_text("413000, 1800000, 950.0\n")
    s = get_session()
    s.model_type = "SRH-2D"
    s.data = object()
    s.project_file = str(tmp_path / "test.srhhydro")

    result = run_calibration(
        [{"type": "manning_n", "material_id": 1, "initial": 0.04,
          "min": 0.02, "max": 0.08}],
        str(csv_file),
        method="unknown_algo",
    )
    assert result["status"] == "error"
    assert "unknown method" in result["message"].lower()
