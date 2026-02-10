from pathlib import Path

import pyHMT2D


def test_monte_carlo_generates_samples_without_solver(tmp_path: Path):
    payload = {
        "srh_control_file": "tests/data/Muncie.srhhydro",
        "n_samples": 2,
        "output_dir": str(tmp_path / "mc_runs"),
        "manning_n_ranges": [{"material_id": 1, "min": 0.02, "max": 0.03}],
        "inlet_q_ranges": [{"bc_id": 1, "min": 50.0, "max": 60.0}],
        "run_simulations": False,
        "seed": 99,
    }

    res = pyHMT2D.Agentic.run_monte_carlo_srh(payload)
    data = res.to_dict()

    assert data["status"] == "succeeded"
    assert data["metrics"]["n_samples"] == 2

    out_dir = Path(payload["output_dir"])
    assert (out_dir / "sample_0000").exists()
    assert (out_dir / "sample_0001").exists()
    assert (out_dir / "monte_carlo_samples.json").exists()
    assert (out_dir / "monte_carlo_samples.csv").exists()
