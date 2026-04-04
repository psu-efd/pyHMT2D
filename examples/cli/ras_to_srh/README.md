## Readme

This directory contains an example for converting a HEC-RAS 2D case to SRH-2D using `hmt-cli`.

Make sure you are in a Windows **PowerShell** terminal with the right Python virtual environment activated. Then, in this example's directory, run:

```powershell
cd Muncie
hmt-cli ras_to_srh --args '{"ras_hdf_file": "Muncie2D.p01.hdf", "srh_case_name": "srh_Muncie"}'
```

The `ras_to_srh` tool takes two required arguments and one optional argument:

- **ras_hdf_file** *(required)*: HEC-RAS 2D result file in HDF format (e.g., `Muncie2D.p01.hdf`).
  The tool derives the project (`.prj`) file and plan ID automatically from this path.
- **srh_case_name** *(required)*: SRH-2D case name prefix (e.g., `srh_Muncie`), which causes the
  command to generate `srh_Muncie.srhgeom` and `srh_Muncie.srhmat`. The `.srhhydro` file needs
  to be created separately.
- **terrain_tif_file** *(optional)*: Terrain file in GeoTIFF format. If omitted, the terrain is
  auto-detected from the geometry HDF (`Terrain Filename` attribute) and the `.rasmap` file.
  Only provide this if auto-detection fails.

> **cmd.exe alternative** (escaped double quotes):
> ```cmd
> hmt-cli ras_to_srh --args "{\"ras_hdf_file\": \"Muncie2D.p01.hdf\", \"srh_case_name\": \"srh_Muncie\"}"
> ```

For a list of all available tools:

```powershell
hmt-cli --help
```
