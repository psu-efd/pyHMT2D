"""
pyHMT2D's CLI: hmt_calibrate

This command line interface calibrates a 2D model. The example command syntax:

$ hmt_calibrate calibration.json

Here, the command "hmt_calibrate" takes one argument:
    - calibration specification JSON file, e.g., "calibration.json"

You can also type the following for more information:

$ hmt_calibrate -v
$ hmt_calibrate -h

"""

import argparse

import pyHMT2D

from ..__about__ import get_pyHMT2D_version_info

def hmt_calibrate(argv=None):
    # Parse command line arguments.
    parser = get_calibrate_parser()
    args = parser.parse_args(argv)

    #create the Calibrator object
    my_calibrator = pyHMT2D.Calibration.Calibrator(args.json_file)

    my_calibrator.calibrate()

    print("Done!")

def get_calibrate_parser():
    parser = argparse.ArgumentParser(
        description=("Calibrate a hydraulic model. The calibration specification file in JSON format needs to be provided as an argument.")
    )

    parser.add_argument("json_file", type=str, help="Calibration specification file in JSON format")

    parser.add_argument(
        "--version",
        "-v",
        action="version",
        version=get_pyHMT2D_version_info(),
        help="Print pyHMT2D version information",
    )

    return parser
