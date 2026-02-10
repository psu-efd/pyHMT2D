import json
from pathlib import Path

import pyHMT2D


def test_validate_job_spec_success():
    payload = {
        "request_type": "result_summary",
        "root": "tests/data",
        "suffixes": [".hdf"],
        "max_files": 10,
    }
    res = pyHMT2D.Agentic.validate_job_spec(payload)
    data = res.to_dict()

    assert data["status"] == "succeeded"
    assert "normalized_request" in data["metrics"]


def test_validate_job_spec_failure():
    payload = {
        "request_type": "convert",
        "ras_hdf_file": "tests/data/Muncie2D.p01.hdf",
        # missing terrain_tif_file and srh_case_name
    }
    res = pyHMT2D.Agentic.validate_job_spec(payload)
    data = res.to_dict()

    assert data["status"] == "failed"
    assert data["error_code"] == "VALIDATION_ERROR"


def test_summarize_artifacts_creates_manifest(tmp_path: Path):
    (tmp_path / "a.json").write_text(json.dumps({"a": 1}), encoding="utf-8")
    (tmp_path / "b.log").write_text("log", encoding="utf-8")

    payload = {
        "root": str(tmp_path),
        "suffixes": [".json", ".log"],
        "max_files": 20,
    }

    res = pyHMT2D.Agentic.summarize_artifacts(payload)
    data = res.to_dict()

    assert data["status"] == "succeeded"
    assert data["metrics"]["file_count"] >= 2
    assert any(Path(a["path"]).name == "run_manifest.json" for a in data["artifacts"])
