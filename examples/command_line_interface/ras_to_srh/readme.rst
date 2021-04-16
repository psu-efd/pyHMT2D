Readme
-----------

This directory contains examples for the "hmt-ras-to-srh" command of *pyHMT2D*.

Make sure you are in a Windows terminal with the right Python virtual environment activated. Then, in this example's
directory, do:

.. code-block:: bash

   $ cd Munice
   $ hmt-ras-to-srh Muncie2D.p01.hdf Terrain/TerrainMuncie_composite.tif srh_Muncie

Here, the command "htm-ras-to-srh" take three arguments:

- HEC-RAS 2D result file in HDF format, e.g., "Muncie2D.p01.hdf"
- Terrain file in GeoTiff format, e.g., "Terrain/TerrainMuncie_composite.tif"
- SRH-2D case name, e.g., "srh_Muncie", which will direct the command to generate "srh_Muncie.srhgeom" and "srh_Muncie.srhmat" files. The ".srhhydro" file needs to be created separately.

You can also type the following for more information:

.. code-block:: bash

   $ htm-ras-to-srh -v
   $ htm-ras-to-srh -h
