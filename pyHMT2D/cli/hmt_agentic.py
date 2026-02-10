"""pyHMT2D agent-facing CLI for typed workflow requests."""

import argparse
import json
from pathlib import Path

import pyHMT2D

from ..__about__ import get_pyHMT2D_version_info


def hmt_agentic(argv=None):
    parser = get_agentic_parser()
    args = parser.parse_args(argv)

    with open(args.request_json, "r", encoding="utf-8") as f:
        payload = json.load(f)

    if args.tool == "validate":
        res = pyHMT2D.Agentic.validate_job_spec(payload)
    elif args.tool == "calibrate":
        res = pyHMT2D.Agentic.run_calibration(payload)
    elif args.tool == "convert":
        res = pyHMT2D.Agentic.convert_ras_to_srh(payload)
    elif args.tool == "summarize":
        res = pyHMT2D.Agentic.summarize_artifacts(payload)
    elif args.tool == "run-model":
        res = pyHMT2D.Agentic.run_single_simulation(payload)
    elif args.tool == "modify-srh":
        res = pyHMT2D.Agentic.modify_srh_case_parameters(payload)
    elif args.tool == "monte-carlo-srh":
        res = pyHMT2D.Agentic.run_monte_carlo_srh(payload)
    else:
        raise ValueError(f"Unsupported tool: {args.tool}")

    result_dict = res.to_dict()
    print(json.dumps(result_dict, indent=2))

    if args.output_json:
        output = Path(args.output_json)
        with open(output, "w", encoding="utf-8") as f:
            json.dump(result_dict, f, indent=2)


def get_agentic_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Run typed pyHMT2D agent-facing tools from JSON request files."
        )
    )

    parser.add_argument(
        "tool",
        choices=["validate", "calibrate", "convert", "summarize", "run-model", "modify-srh", "monte-carlo-srh"],
        help="Tool to execute",
    )
    parser.add_argument("request_json", type=str, help="Path to JSON request payload")
    parser.add_argument("--output-json", type=str, default=None, help="Optional output JSON path")

    parser.add_argument(
        "--version",
        "-v",
        action="version",
        version=get_pyHMT2D_version_info(),
        help="Print pyHMT2D version information",
    )

    return parser
