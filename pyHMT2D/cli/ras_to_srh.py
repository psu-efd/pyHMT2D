import argparse

import pyHMT2D

from ..__about__ import get_pyHMT2D_version_info

def ras_to_srh(argv=None):
    # Parse command line arguments.
    parser = get_ras_to_srh_parser()
    args = parser.parse_args(argv)

    #build the RAS_to_SRH_Converter object
    my_ras_to_srh_converter = pyHMT2D.Misc.RAS_to_SRH_Converter(args.ras_infile,
                                                                args.ras_terrainfile,
                                                                args.srh_casename)

    #convert to SRH-2D (mesh and material) and save
    my_ras_to_srh_converter.convert_to_SRH()

    print("Done!")

def get_ras_to_srh_parser():
    parser = argparse.ArgumentParser(
        description=("Convert HEC-RAS 2D mesh and Manning's n to SRH-2D. It will create two files: .srhgeom and .srhmat.")
    )

    parser.add_argument("ras_infile", type=str, help="HEC-RAS 2D result hdf file, e.g., case.p01.hdf")

    parser.add_argument("ras_terrainfile", type=str, help="HEC-RAS 2D terrain file, e.g., Terrain/terrain.tif")

    parser.add_argument("srh_casename", type=str, help="SRH-2D case name. The result files will be"
                                                       " case_name.srhgeom and case_name.srhmat")

    parser.add_argument(
        "--version",
        "-v",
        action="version",
        version=get_pyHMT2D_version_info(),
        help="Print pyHMT2D version information",
    )

    return parser
