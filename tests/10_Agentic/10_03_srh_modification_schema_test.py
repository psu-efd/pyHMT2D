import pyHMT2D


def test_modify_srh_request_validation_success():
    payload = {
        "request_type": "modify_srh_case",
        "srh_control_file": "tests/data/Muncie.srhhydro",
        "material_ids": [1, 2],
        "manning_n_values": [0.03, 0.04],
        "inlet_bc_ids": [1],
        "inlet_q_values": [500.0],
    }

    res = pyHMT2D.Agentic.validate_job_spec(payload)
    data = res.to_dict()

    assert data["status"] == "succeeded"
    assert data["metrics"]["normalized_request"]["material_ids"] == [1, 2]


def test_monte_carlo_request_validation_success(tmp_path):
    payload = {
        "request_type": "monte_carlo_srh",
        "srh_control_file": "tests/data/Muncie.srhhydro",
        "n_samples": 3,
        "output_dir": str(tmp_path / "mc"),
        "manning_n_ranges": [{"material_id": 1, "min": 0.02, "max": 0.05}],
        "inlet_q_ranges": [{"bc_id": 1, "min": 100.0, "max": 300.0}],
        "run_simulations": False,
    }

    res = pyHMT2D.Agentic.validate_job_spec(payload)
    data = res.to_dict()

    assert data["status"] == "succeeded"
    assert data["metrics"]["normalized_request"]["n_samples"] == 3
