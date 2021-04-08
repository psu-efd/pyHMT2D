# Comparison between SRH-2D and HEC-RAS 2D

This directory contains cases for the comparison between SRH-2D and HEC-RAS. The workflow is typicallly as follows:

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
- Process RAS 2D data:
    - convert RAS 2D mesh and material (Manning's n) with the RAS_to_SRH_Converter class. This will generate "srhgeom" 
      and "srhmat" files. The srhhydro file needs to be created manually (see SRH-2D and SMS examples for reference)
    - convert RAS 2D result to VTK (for further processing and comparsion with SRH-2D)
-     


## Simulation cases:
- Backwater curves
- Muncie

    