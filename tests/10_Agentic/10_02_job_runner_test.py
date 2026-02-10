from pathlib import Path

import pyHMT2D


def test_run_job_persists_manifest(tmp_path: Path):
    payload = {
        "request_type": "result_summary",
        "root": "tests/data",
        "suffixes": [".hdf"],
        "max_files": 3,
    }

    handle = pyHMT2D.Agentic.run_job(
        job_dir=str(tmp_path),
        job_type="validate",
        fn=pyHMT2D.Agentic.validate_job_spec,
        payload=payload,
    )

    assert handle.manifest_path.exists()
    content = handle.manifest_path.read_text(encoding="utf-8")
    assert '"job_type": "validate"' in content
