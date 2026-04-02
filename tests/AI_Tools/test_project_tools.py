# -*- coding: utf-8 -*-
"""Tests for AI_Tools project_tools (no solver required)."""

import os
import pytest

from pyHMT2D.AI_Tools.state import reset_session, get_session
from pyHMT2D.AI_Tools.tools.project_tools import (
    list_model_files,
    get_project_info,
    get_materials,
    get_boundary_conditions,
    get_result_variables,
)


@pytest.fixture(autouse=True)
def clean_session():
    """Reset the global session before each test."""
    reset_session()
    yield
    reset_session()


# ── list_model_files ──────────────────────────────────────────────────────────

def test_list_model_files_missing_dir():
    result = list_model_files("/nonexistent/path/xyz")
    assert result["status"] == "error"
    assert "not found" in result["message"].lower()


def test_list_model_files_empty_dir(tmp_path):
    result = list_model_files(str(tmp_path))
    assert result["status"] == "ok"
    assert result["data"] == []


def test_list_model_files_detects_srhhydro(tmp_path):
    (tmp_path / "Muncie.srhhydro").write_text("SRH2D test")
    result = list_model_files(str(tmp_path))
    assert result["status"] == "ok"
    assert len(result["data"]) == 1
    assert result["data"][0]["model_type"] == "SRH-2D"
    assert result["data"][0]["case_name"] == "Muncie"


def test_list_model_files_detects_hdf(tmp_path):
    (tmp_path / "Muncie2D.p01.hdf").write_bytes(b"\x89HDF")
    result = list_model_files(str(tmp_path))
    assert result["status"] == "ok"
    assert any(f["model_type"] == "HEC-RAS" for f in result["data"])


# ── get_project_info ──────────────────────────────────────────────────────────

def test_get_project_info_missing_file():
    result = get_project_info("/no/such/file.srhhydro")
    assert result["status"] == "error"
    assert "not found" in result["message"].lower()


def test_get_project_info_unknown_extension(tmp_path):
    f = tmp_path / "model.xyz"
    f.write_text("data")
    result = get_project_info(str(f))
    assert result["status"] == "error"
    assert "cannot determine model type" in result["message"].lower()


# ── session state after get_project_info ─────────────────────────────────────

def test_require_open_before_get_materials():
    result = get_materials()
    assert result["status"] == "error"
    assert "no project" in result["message"].lower()


def test_require_open_before_get_boundary_conditions():
    result = get_boundary_conditions()
    assert result["status"] == "error"
    assert "no project" in result["message"].lower()


# ── get_result_variables ──────────────────────────────────────────────────────

def test_get_result_variables_missing_file():
    result = get_result_variables("/no/such/file.h5")
    assert result["status"] == "error"
    assert "not found" in result["message"].lower()
