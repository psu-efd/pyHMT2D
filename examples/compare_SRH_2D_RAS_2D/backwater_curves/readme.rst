Readme
==============================================

To perform the comparison, the order of steps is as follows:

- Run "create_channel_terrain.py": it will create the constant slope channel terrain in GeoTiff format, which can be used in HEC-RAS for case building.
- In folder "HEC-RAS-2D", build the simple 2D case with the created terrain (the case has been provided). Run the case to get result.
- Run "process_RAS_2D_Data.py": it will process HEC-RAS 2D result and save to VTK. It will also convert HEC-RAS 2D mesh and Manning's n to SRH-2D.
- In folder "SRH-2D", run the SRH-2D case to get result.
- Run "process_SRH_2D_Data.py": it will process SRH-2D result and save to VTK.
- Run "compare_SRH_2D_HEC_RAS_2D.py": it will read the two VTK result files from both SRH-2D and HEC-RAS 2D, sample \
  the water surface profile, run a simple Backwater-1D solver, and finally plot all profiles in one figure. \
  It also calculates the difference in results between SRH-2D and HEC-RAS 2D and save to VTk files, which can be \
  loaded in ParaView for inspection.

    .. figure:: backwater_1D_comparison.png
        :width: 400px
        :align: center
        :height: 320px
        :alt: backwater curve
        :figclass: align-center


    .. figure:: backwater_diff_paraview.png
        :width: 400px
        :align: center
        :height: 320px
        :alt: backwater curve
        :figclass: align-center


