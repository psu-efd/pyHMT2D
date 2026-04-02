## Readme

This directory contains an example for exporting only the SRH-2D **mesh** to VTK using `hmt-cli`.

Make sure you are in a Windows **PowerShell** terminal with the correct Python virtual environment activated. Then, in this example's directory, run:

```powershell
hmt-cli srh_to_vtk --args '{"srhhydro_file": "Muncie.srhhydro", "output_file": "my_case_mesh.vtk", "mesh_only": true}'
```

Arguments:

- **srhhydro_file**: SRH-2D case file in `.srhhydro` format (e.g., `Muncie.srhhydro`).
- **output_file**: Output VTK file name; must include the `.vtk` extension (e.g., `my_case_mesh.vtk`).
- **mesh_only**: Set to `true` to export only the mesh with no result data.

> **cmd.exe alternative** (escaped double quotes):
> ```cmd
> hmt-cli srh_to_vtk --args "{\"srhhydro_file\": \"Muncie.srhhydro\", \"output_file\": \"my_case_mesh.vtk\", \"mesh_only\": true}"
> ```
