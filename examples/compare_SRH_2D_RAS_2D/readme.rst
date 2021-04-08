Comparison between SRH-2D and HEC-RAS 2D
========================================

General steps:
##############

This directory contains cases for the comparison between SRH-2D and HEC-RAS using the same mesh. Indeed, HEC-RAS is \
used as a mesh generator for SRH-2D.

The workflow is typicallly as follows:

- Prepare data files such as bathymetry for HEC-RAS, for example 
    - use the Terrain class to create terrain GeoTiff file
- In HEC-RAS, create a case:
    - load terrain data
    - create geometry with only one 2D flow area and generate the mesh
    - create boundary condition
    - create Manning's n layer if Manning's n is spatially varying
    - create unsteady flow
    - create simulation plan
    - run simulation
- Process HEC-RAS 2D data:
    - convert RAS 2D mesh and material (Manning's n) with the RAS_to_SRH_Converter class. This will generate "srhgeom" \
      and "srhmat" files. The srhhydro file needs to be created manually (see SRH-2D and SMS examples for reference)
    - convert RAS 2D result to VTK (for further processing and comparsion with SRH-2D)
- In SRH-2D,
    - make necessary modifications to the SRH-2D case files (srhhydro, srhgeom, srhmat, etc.)
    - run SRH-2D preprocessing and solver
- Process SRH-2D data:
    - convert SRH-2D result to VTk
- Compare and plot results
    - read in the VTK result files from both SRH-2D and HEC-RAS 2D
    - extract profiles, probe on points, calculate differences, etc.
    - plot


Example cases:
##############

- Backwater curve
- Muncie

    