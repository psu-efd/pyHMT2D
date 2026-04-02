## Readme

--------------------------------
**Note**: This example requires SRH-2D results in HDF format. You need to run the SRH-2D simulation first to generate the HDF file. For example, run the following command in the SRH-2D directory (assuming you are using Aquaveo's SMS 13.4):
```cmd
"C:\Program Files\SMS 13.4 64-bit\python\Lib\site-packages\srh2d_exe\SRH_Pre_Console.exe" 3 Muncie.srhhydro
"C:\Program Files\SMS 13.4 64-bit\python\Lib\site-packages\srh2d_exe\SRH-2D_Console.exe" Muncie.DAT
```
This will generate the HDF file `Muncie_XMDFC.h5` in the same directory.

--------------------------------

This directory contains an example for converting SRH-2D **results** (solution data from an HDF file) to VTK using `hmt-cli`.

Make sure you are in a Windows **PowerShell** terminal with the correct Python virtual environment activated. Then, in this example's directory, run:

```powershell
hmt-cli srh_to_vtk --args '{"srhhydro_file": "Muncie.srhhydro", "output_file": "Muncie_XMDFC.h5"}'
```

Arguments:

- **srhhydro_file**: SRH-2D case file in `.srhhydro` format (e.g., `Muncie.srhhydro`).
- **output_file**: SRH-2D result file in HDF format (e.g., `Muncie_XMDFC.h5`). The command reads this file and exports the solution to VTK.

Optional fields:

- **nodal** (`true`/`false`): use nodal (node) data instead of cell-centered (default: `false`).
- **all_timesteps** (`true`/`false`): convert all time steps to VTK; if not set, only the last time step is converted (default: `false`).

Example with options:

```powershell
hmt-cli srh_to_vtk --args '{"srhhydro_file": "Muncie.srhhydro", "output_file": "Muncie_XMDFC.h5", "nodal": true, "all_timesteps": true}'
```

> **cmd.exe alternative** (escaped double quotes):
> ```cmd
> hmt-cli srh_to_vtk --args "{\"srhhydro_file\": \"Muncie.srhhydro\", \"output_file\": \"Muncie_XMDFC.h5\"}"
> ```

To export only the mesh (no result data), see the `srh_mesh_to_vtk` example directory.
