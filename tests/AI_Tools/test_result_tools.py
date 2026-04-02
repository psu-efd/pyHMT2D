# -*- coding: utf-8 -*-
"""Tests for AI_Tools result_tools (no solver required)."""

import pytest

from pyHMT2D.AI_Tools.state import reset_session
from pyHMT2D.AI_Tools.tools.result_tools import (
    read_results,
    get_value_at_point,
    get_result_statistics,
    get_flood_extent,
    get_cross_section_profile,
)


@pytest.fixture(autouse=True)
def clean_session():
    reset_session()
    yield
    reset_session()


def test_read_results_missing_file():
    from pyHMT2D.AI_Tools.state import get_session
    s = get_session()
    s.model_type = "SRH-2D"
    s.data = object()

    result = read_results("/no/such/file.h5")
    assert result["status"] == "error"
    assert "not found" in result["message"].lower()


def test_get_value_at_point_no_results():
    result = get_value_at_point(0.0, 0.0, "Water_Elev_m")
    assert result["status"] == "error"
    assert "no project" in result["message"].lower() or "results" in result["message"].lower()


def test_get_result_statistics_no_results():
    result = get_result_statistics("Water_Depth_m")
    assert result["status"] == "error"


def test_get_flood_extent_no_results():
    result = get_flood_extent(depth_threshold=0.01)
    assert result["status"] == "error"


def test_get_cross_section_profile_no_results():
    result = get_cross_section_profile(0, 0, 100, 100, "Water_Elev_m")
    assert result["status"] == "error"


def test_get_cross_section_profile_bad_n_points():
    from pyHMT2D.AI_Tools.state import get_session
    s = get_session()
    s.model_type = "SRH-2D"
    s.data = object()
    s.result_loaded = True

    result = get_cross_section_profile(0, 0, 100, 100, "Water_Elev_m", n_points=1)
    assert result["status"] == "error"
    assert "at least 2" in result["message"]
