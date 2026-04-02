# -*- coding: utf-8 -*-
"""Tests for AI_Tools parameter_tools (no solver required)."""

import pytest

from pyHMT2D.AI_Tools.state import reset_session
from pyHMT2D.AI_Tools.tools.parameter_tools import (
    set_manning_n,
    set_inlet_flow,
    set_exit_wse,
    save_modified_inputs,
)


@pytest.fixture(autouse=True)
def clean_session():
    reset_session()
    yield
    reset_session()


# ── set_manning_n ─────────────────────────────────────────────────────────────

def test_set_manning_n_no_project():
    result = set_manning_n([1], [0.04])
    assert result["status"] == "error"
    assert "no project" in result["message"].lower()


def test_set_manning_n_mismatched_lengths():
    # Inject a fake open session
    from pyHMT2D.AI_Tools.state import get_session
    s = get_session()
    s.model_type = "SRH-2D"
    s.data = object()  # placeholder

    result = set_manning_n([1, 2], [0.04])  # 2 IDs, 1 value
    assert result["status"] == "error"
    assert "same length" in result["message"]


def test_set_manning_n_out_of_range():
    from pyHMT2D.AI_Tools.state import get_session
    s = get_session()
    s.model_type = "SRH-2D"
    s.data = object()

    result = set_manning_n([1], [5.0])  # 5.0 > 1.0 physical limit
    assert result["status"] == "error"
    assert "plausible range" in result["message"]


# ── set_inlet_flow ────────────────────────────────────────────────────────────

def test_set_inlet_flow_no_project():
    result = set_inlet_flow([1], [100.0])
    assert result["status"] == "error"


def test_set_inlet_flow_hecras_rejected():
    from pyHMT2D.AI_Tools.state import get_session
    s = get_session()
    s.model_type = "HEC-RAS"
    s.data = object()

    result = set_inlet_flow([1], [100.0])
    assert result["status"] == "error"
    assert "srh-2d" in result["message"].lower()


def test_set_inlet_flow_negative_value():
    from pyHMT2D.AI_Tools.state import get_session
    s = get_session()
    s.model_type = "SRH-2D"
    s.data = object()

    result = set_inlet_flow([1], [-50.0])
    assert result["status"] == "error"
    assert "non-negative" in result["message"]


# ── set_exit_wse ──────────────────────────────────────────────────────────────

def test_set_exit_wse_hecras_rejected():
    from pyHMT2D.AI_Tools.state import get_session
    s = get_session()
    s.model_type = "HEC-RAS"
    s.data = object()

    result = set_exit_wse([1], [950.0])
    assert result["status"] == "error"


# ── save_modified_inputs ──────────────────────────────────────────────────────

def test_save_modified_inputs_no_project():
    result = save_modified_inputs()
    assert result["status"] == "error"
