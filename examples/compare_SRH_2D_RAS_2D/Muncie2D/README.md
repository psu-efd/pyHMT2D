# Readme

This is a more realistic case using the terrain data in the "Muncie" example that comes with HEC-RAS. The domain is a simple 2D flow area with one inlet (fixed discharge) and one outlet (fixed water surface elevation). Other things such as hydraulic structures and 1D channels can exist in the HEC-RAS case. However, **pyHMT2D**'s RAS_2D_Data class currently does not support these. Thus, if you have these, you may have to add them in the SRH-2D case manually.

To perform the comparison, the order of steps is as follows:

- In folder "`HEC-RAS`", load the case in HEC-RAS (version 6.6) and run the case to get result.
- Run "`python process_RAS_2D_Data.py`": it will process HEC-RAS 2D result and save to VTK. It will also convert HEC-RAS 2D mesh and Manning's $n$ to SRH-2D, and result in two SRH-2D files: `Muncie.srhgeom` and `Muncie.srhmat`.
- Copy the `Muncie.srhgeom` and `Muncie.srhmat` files to the `SRH-2D` folder. SRH-2D also requires the `Muncie.srhhydro` file, which has already been provided in the `SRH-2D` folder.
- In folder "`SRH-2D`", run the SRH-2D case to get result by running the following commands (assuming you are using Aquaveo's SMS 13.4):
  
```bash
# run SRH-2D pre-processor to generate the DAT file
"C:\Program Files\SMS 13.4 64-bit\python\Lib\site-packages\srh2d_exe\SRH_Pre_Console.exe" 3 Muncie.srhhydro

# run SRH-2D to generate the HDF file
"C:\Program Files\SMS 13.4 64-bit\python\Lib\site-packages\srh2d_exe\SRH-2D_Console.exe" Muncie.DAT
```


- Run "`python process_SRH_2D_Data.py`": it will process SRH-2D result and save to VTK.
- Run "`python compare_SRH_2D_HEC_RAS_2D.py`": it will read the two VTK result files from both SRH-2D and HEC-RAS 2D, calculates the difference in results between SRH-2D and HEC-RAS 2D, and save to VTK files, which can be loaded in ParaView for inspection.

Comparison of water depth simulated by HEC-RAS and SRH-2D:

![Muncie comparison between SRH-2D and HEC-RAS](Muncie_Paraview_compare_water_depth.png)

Difference between HEC-RAS and SRH-2D in bed elevation, water depth, velocity, and WSE:

![Difference between SRH-2D and HEC-RAS](Muncie_differences_paraview.png)
