"""
pyHMT2D's CLI: hmt-srh-mesh-to-vtk

This command line interface converts a SRH-2D mesh to VTK.

The example command syntax:

$ hmt-srh-mesh-to-vtk case.srhhydro case_mesh.vtk

Here, the command "hmt-srh-mesh-to-vtk" takes two arguments:
    - SRH-2D's srhhydro file name, e.g., "case.srhhydro"
    - Output VTK file name, e.g., "case_mesh.vtk"

You can also type the following for more information:

$ hmt-srh-mesh-to-vtk -v
$ hmt-srh-mesh-to-vtk -h

"""

import argparse

import pyHMT2D

from ..__about__ import get_pyHMT2D_version_info

def hmt_srh_mesh_to_vtk(argv=None):
    # Parse command line arguments.
    parser = get_srh_mesh_to_vtk_parser()
    args = parser.parse_args(argv)

    # read the SRH-2D's case files (srhhydro, srhgeom, and srhmat)
    my_srh_2d_data = pyHMT2D.SRH_2D.SRH_2D_Data(args.srhhydro_file)



    #my_srh_2d_data.ManningN_cell

    #my_srh_2d_data.outputVTK(vtkFileName, resultVarNames, resultData, bNodal):

    # output the 3D VTK file (mesh not flat)
    my_srh_2d_data.srhgeom_obj.output_2d_mesh_to_vtk(args.vtk_file, bFlat=False)

    print("Done!")

def get_srh_mesh_to_vtk_parser():
    parser = argparse.ArgumentParser(
        description=("Converts a SRH-2D mesh to VTK.")
    )

    parser.add_argument("srhhydro_file", type=str, help="SRH-2D case file in srhhydro format")

    parser.add_argument("vtk_file", type=str, help="Output VTK file name; must include .vtk")

    parser.add_argument(
        "--version",
        "-v",
        action="version",
        version=get_pyHMT2D_version_info(),
        help="Print pyHMT2D version information",
    )

    return parser
