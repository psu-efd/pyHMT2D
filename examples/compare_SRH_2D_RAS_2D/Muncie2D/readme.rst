Readme
==============================================

This is a more realistic case using the terrain data in the "Muncie" example that comes with HEC-RAS. The domain is a simple 2D flow area with one inlet (fixed discharge) and one outlet (fixed water surface elevation). Other things such as hydraulic structures and 1D channels can exist in the HEC-RAS case. However, *pyHMT2D*'s RAS_2D_Data class currently does not support these. Thus, if you have these, you may have to add them in the SRH-2D case manually.

To perform the comparison, the order of steps is as follows:

- In folder "HEC-RAS", build the simple 2D case with the created terrain (the case has been provided). Run the case to get result.
- Run "process_RAS_2D_Data.py": it will process HEC-RAS 2D result and save to VTK. It will also convert HEC-RAS 2D mesh and Manning's n to SRH-2D.
- In folder "SRH-2D", run the SRH-2D case to get result.
- Run "process_SRH_2D_Data.py": it will process SRH-2D result and save to VTK.
- Run "compare_SRH_2D_HEC_RAS_2D.py": it will read the two VTK result files from both SRH-2D and HEC-RAS 2D, calculates the difference in results between SRH-2D and HEC-RAS 2D, and save to VTk files, which can be loaded in ParaView for inspection.

Comparison of water depth simulated by HEC-RAS and SRH-2D:

    .. figure:: Muncie_Paraview_compare_water_depth.png
        :width: 400px
        :align: center
        :height: 320px
        :alt: Muncie comparison between SRH-2D and HEC-RAS
        :figclass: align-center


Difference between HEC-RAS and SRH-2D in bed elevation, water depth, velocity, and WSE:

    .. figure:: Muncie_differences_paraview.png
        :width: 400px
        :align: center
        :height: 320px
        :alt: Difference between SRH-2D and HEC-RAS
        :figclass: align-center

