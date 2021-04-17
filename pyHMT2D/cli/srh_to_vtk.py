"""
pyHMT2D's CLI: hmt-srh-to-vtk

This command line interface converts a SRH-2D result to VTK.

The example command syntax:

$ hmt-srh-to-vtk case.srhhydro case_XMDFC.h5

Here, the command "hmt-srh-mesh-to-vtk" takes two arguments:
    - SRH-2D's srhhydro file name, e.g., "case.srhhydro"
    - Output VTK file name, e.g., "case_mesh.vtk"

You can also type the following for more information:

$ hmt-srh-to-vtk -v
$ hmt-srh-to-vtk -h

"""

import argparse

import pyHMT2D

from ..__about__ import get_pyHMT2D_version_info

def hmt_srh_to_vtk(argv=None):
    # Parse command line arguments.
    parser = get_srh_to_vtk_parser()
    args = parser.parse_args(argv)

    # read the SRH-2D's case files (srhhydro, srhgeom, and srhmat)
    my_srh_2d_data = pyHMT2D.SRH_2D.SRH_2D_Data(args.srhhydro_file)

    #read SRH-2D result in XMDF format (*.h5)
    my_srh_2d_data.readSRHXMDFFile(args.hdf_file, bNodal=False)

    #export to VTK
    my_srh_2d_data.outputXMDFDataToVTK(bNodal=False, lastTimeStep=True, dir='')

    print("Done!")

def get_srh_to_vtk_parser():
    parser = argparse.ArgumentParser(
        description=("Converts a SRH-2D result to VTK.")
    )

    parser.add_argument("srhhydro_file", type=str, help="SRH-2D case file in srhhydro format")

    parser.add_argument("hdf_file", type=str, help="SRH-2D result file in HDF format; must include .h5")

    parser.add_argument(
        "--version",
        "-v",
        action="version",
        version=get_pyHMT2D_version_info(),
        help="Print pyHMT2D version information",
    )

    return parser
