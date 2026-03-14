"""
pyHMT2D's CLI: hmt-srh-to-vtk

This command line interface converts a SRH-2D result to VTK, or exports only the mesh.

Examples:

    # Result to VTK (reads HDF result, exports solution to VTK)
    $ hmt-srh-to-vtk case.srhhydro case_XMDFC.h5

    # Mesh only (no HDF; export mesh to given VTK file)
    $ hmt-srh-to-vtk --mesh-only case.srhhydro case_mesh.vtk

Arguments:
    - srhhydro_file: SRH-2D case file in srhhydro format
    - hdf_or_vtk_file: HDF result file (.h5) for result conversion; or output VTK file when --mesh-only

Optional:
    - --mesh-only: export only the SRH-2D mesh to VTK (no result data)
    - --nodal: use nodal data instead of cell-centered (default: False; ignored for --mesh-only)
    - --all-timesteps: convert all time steps; default is last only (ignored for --mesh-only)

    $ hmt-srh-to-vtk -v
    $ hmt-srh-to-vtk -h
"""

import argparse

import pyHMT2D

from ..__about__ import get_pyHMT2D_version_info


def hmt_srh_to_vtk(argv=None):
    parser = get_srh_to_vtk_parser()
    args = parser.parse_args(argv)

    my_srh_2d_data = pyHMT2D.SRH_2D.SRH_2D_Data(args.srhhydro_file)

    if args.mesh_only:
        # Export only the mesh to VTK
        my_srh_2d_data.srhgeom_obj.output_2d_mesh_to_vtk(args.hdf_or_vtk_file, bFlat=False)
    else:
        # Read HDF result and export to VTK
        my_srh_2d_data.readSRHXMDFFile(args.hdf_or_vtk_file, bNodal=args.nodal)
        my_srh_2d_data.outputXMDFDataToVTK(
            bNodal=args.nodal,
            lastTimeStep=not args.all_timesteps,
            dir='',
        )

    print("Done!")

def get_srh_to_vtk_parser():
    parser = argparse.ArgumentParser(
        description="Converts a SRH-2D result to VTK, or exports only the mesh (--mesh-only).",
    )

    parser.add_argument("srhhydro_file", type=str, help="SRH-2D case file in srhhydro format")

    parser.add_argument(
        "hdf_or_vtk_file",
        type=str,
        help="HDF result file (.h5) for result-to-VTK; or output VTK file when using --mesh-only",
    )

    parser.add_argument(
        "--mesh-only",
        action="store_true",
        default=False,
        help="Export only the SRH-2D mesh to VTK (second argument is output .vtk file)",
    )

    parser.add_argument(
        "--nodal",
        action="store_true",
        default=False,
        help="Use nodal (node) data instead of cell-centered (default: False; ignored if --mesh-only)",
    )

    parser.add_argument(
        "--all-timesteps",
        action="store_true",
        default=False,
        help="Convert all time steps to VTK; if not set, only last time step (default: last only; ignored if --mesh-only)",
    )

    parser.add_argument(
        "--version",
        "-v",
        action="version",
        version=get_pyHMT2D_version_info(),
        help="Print pyHMT2D version information",
    )

    return parser
