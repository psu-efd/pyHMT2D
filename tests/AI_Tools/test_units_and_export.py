# -*- coding: utf-8 -*-
"""Tests for unit resolution and VTK export bookkeeping (no solver required)."""

import os
from types import SimpleNamespace

import pytest

from pyHMT2D.AI_Tools.state import get_session, reset_session
from pyHMT2D.AI_Tools.tools._units import (
    area_to_km2,
    area_unit,
    discharge_unit,
    is_english,
    length_unit,
)
from pyHMT2D.AI_Tools.tools.export_tools import _vtk_fingerprints, _written_since
from pyHMT2D.AI_Tools.tools.parameter_tools import set_inlet_flow


@pytest.fixture(autouse=True)
def clean_session():
    reset_session()
    yield
    reset_session()


# ── helpers ───────────────────────────────────────────────────────────────────

def _srh_session(output_unit):
    """A stand-in session for an SRH-2D project reporting the given OutputUnit."""
    data = SimpleNamespace(
        get_output_format_unit=lambda: ("XMDFC", output_unit)
    )
    return SimpleNamespace(model_type="SRH-2D", data=data)


def _ras_session(units):
    return SimpleNamespace(model_type="HEC-RAS", data=SimpleNamespace(units=units))


def _write(path, text):
    with open(path, "w") as fh:
        fh.write(text)


# ── unit resolution ───────────────────────────────────────────────────────────

@pytest.mark.parametrize("output_unit", ["EN", "en", "English"])
def test_srh_english_units(output_unit):
    s = _srh_session(output_unit)
    assert is_english(s)
    assert length_unit(s) == "ft"
    assert area_unit(s) == "ft2"
    assert discharge_unit(s) == "ft3/s"


def test_srh_si_units():
    s = _srh_session("SI")
    assert not is_english(s)
    assert length_unit(s) == "m"
    assert area_unit(s) == "m2"
    assert discharge_unit(s) == "m3/s"


def test_hecras_feet_and_meter():
    assert length_unit(_ras_session("Feet")) == "ft"
    assert length_unit(_ras_session("Meter")) == "m"


def test_unknown_and_closed_sessions_fall_back_to_si():
    assert length_unit(_ras_session("Unknown units")) == "m"
    assert length_unit(SimpleNamespace(model_type=None, data=None)) == "m"


def test_units_survive_a_raising_project():
    """A project that cannot report its units must not break the caller."""
    def boom():
        raise RuntimeError("no OutputUnit keyword")

    s = SimpleNamespace(
        model_type="SRH-2D", data=SimpleNamespace(get_output_format_unit=boom)
    )
    assert length_unit(s) == "m"


def test_area_conversion_uses_the_project_unit():
    # 1 acre = 43,560 ft² = 4046.8564224 m²
    acre_ft2 = 43560.0
    assert area_to_km2(acre_ft2, _srh_session("EN")) == pytest.approx(
        4046.8564224 / 1e6
    )
    # SI areas are already m², so the conversion is a plain divide.
    assert area_to_km2(2.5e6, _srh_session("SI")) == pytest.approx(2.5)


# ── export bookkeeping ────────────────────────────────────────────────────────

def test_written_since_ignores_pre_existing_files(tmp_path):
    """The core bug: a later export must not claim earlier files as its own."""
    stale = tmp_path / "case_0006.vtk"
    _write(stale, "old")

    before = _vtk_fingerprints(str(tmp_path))

    fresh = tmp_path / "case_0001.vtk"
    _write(fresh, "new")

    assert _written_since(str(tmp_path), before) == [str(fresh)]


def test_written_since_detects_overwritten_files(tmp_path):
    """A rewrite under the same name still counts as written."""
    target = tmp_path / "case_0001.vtk"
    _write(target, "first")
    before = _vtk_fingerprints(str(tmp_path))

    _write(target, "second pass, different length")
    # Make the change unambiguous regardless of filesystem timestamp resolution.
    os.utime(target, ns=(0, 0))

    assert _written_since(str(tmp_path), before) == [str(target)]


def test_written_since_empty_directory(tmp_path):
    assert _written_since(str(tmp_path), {}) == []


def test_written_since_falls_back_when_nothing_changed(tmp_path):
    """No detectable change must not be reported as 'exported nothing'."""
    existing = tmp_path / "case_0001.vtk"
    _write(existing, "unchanged")
    before = _vtk_fingerprints(str(tmp_path))

    assert _written_since(str(tmp_path), before) == [str(existing)]


def test_fingerprints_skip_non_vtk(tmp_path):
    _write(tmp_path / "case.vtk", "x")
    _write(tmp_path / "notes.txt", "x")
    assert list(_vtk_fingerprints(str(tmp_path))) == [str(tmp_path / "case.vtk")]


# ── set_inlet_flow reporting ──────────────────────────────────────────────────

class _FakeSRHHydro:
    def __init__(self, discharges):
        self._q = discharges

    def get_InletQ(self):
        return self._q


class _FakeSRHData:
    """Minimal SRH_2D_Data stand-in for the inlet-flow reporting path."""

    control_type = "SRHHydro"

    def __init__(self, discharges, output_unit="EN"):
        self.srhhydro_obj = _FakeSRHHydro(discharges)
        self._output_unit = output_unit
        self.applied = None

    def get_output_format_unit(self):
        return ("XMDFC", self._output_unit)

    def modify_InletQ(self, ids, values):
        self.applied = (ids, values)


def test_set_inlet_flow_reports_old_value_and_project_units():
    s = get_session()
    s.model_type = "SRH-2D"
    s.data = _FakeSRHData({2: {"discharge": 20000.0}})

    result = set_inlet_flow([2], [25000.0])

    assert result["status"] == "ok"
    assert result["data"] == [
        {"id": 2, "old_value": 20000.0, "new_value": 25000.0, "units": "ft3/s"}
    ]
    assert s.data.applied == ([2], [25000.0])


def test_set_inlet_flow_units_follow_an_si_project():
    s = get_session()
    s.model_type = "SRH-2D"
    s.data = _FakeSRHData({1: {"discharge": 500.0}}, output_unit="SI")

    result = set_inlet_flow([1], [600.0])
    assert result["data"][0]["units"] == "m3/s"


def test_set_inlet_flow_unknown_bc_reports_nan_old_value():
    s = get_session()
    s.model_type = "SRH-2D"
    s.data = _FakeSRHData({2: {"discharge": 20000.0}})

    result = set_inlet_flow([7], [100.0])
    assert result["status"] == "ok"
    assert result["data"][0]["old_value"] != result["data"][0]["old_value"]  # NaN
