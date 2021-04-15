Readme
-----------

This directory contains examples for command line interface (CLI) of *pyHMT2D*, i.e., directly using commands in
a windows terminal.

Make sure you are in the Python virtual environment that *pyHMT2D* is installed. Then, navigate to the example
directory, for example:

.. code-block:: bash

   $ cd ras_to_srh/Munice
   $ ras_to_srh Muncie2D.p01.hdf Terrain/TerrainMuncie_composite.tif srh_Muncie

Here, the command "ras_to_srh" take three arguments:

- HEC-RAS 2D result file in HDF format, e.g., "Muncie2D.p01.hdf"
- Terrain file in GeoTiff format, e.g., "Terrain/TerrainMuncie_composite.tif"
- SRH-2D case name, e.g., "srh_Muncie", which will direct the command to generate "srh_Muncie.srhgeom" and "srh_Muncie.srhmat" files. The ".srhhydro" file needs to be created separately.

