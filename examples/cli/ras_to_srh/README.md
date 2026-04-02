## Readme

This directory contains an example for converting a HEC-RAS 2D case to SRH-2D using `hmt-cli`.

Make sure you are in a Windows **PowerShell** terminal with the right Python virtual environment activated. Then, in this example's directory, run:

```powershell
cd Muncie
hmt-cli ras_to_srh --args '{"ras_hdf_file": "Muncie2D.p01.hdf", "terrain_tif_file": "Terrain/TerrainMuncie_composite.tif", "srh_case_name": "srh_Muncie"}'
```

The `ras_to_srh` tool takes three arguments:

- **ras_hdf_file**: HEC-RAS 2D result file in HDF format (e.g., `Muncie2D.p01.hdf`)
- **terrain_tif_file**: Terrain file in GeoTIFF format (e.g., `Terrain/TerrainMuncie_composite.tif`)
- **srh_case_name**: SRH-2D case name (e.g., `srh_Muncie`), which causes the command to generate `srh_Muncie.srhgeom` and `srh_Muncie.srhmat`. The `.srhhydro` file needs to be created separately.

> **cmd.exe alternative** (escaped double quotes):
> ```cmd
> hmt-cli ras_to_srh --args "{\"ras_hdf_file\": \"Muncie2D.p01.hdf\", \"terrain_tif_file\": \"Terrain/TerrainMuncie_composite.tif\", \"srh_case_name\": \"srh_Muncie\"}"
> ```

For a list of all available tools:

```powershell
hmt-cli --help
```
