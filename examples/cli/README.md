## Readme

This directory contains examples for the command line interface (CLI) of *pyHMT2D*.
All commands are accessed through the unified `hmt-cli` entry point.

### Shell requirements

`hmt-cli` passes arguments as JSON strings. **Use PowerShell** — single-quoted JSON strings work correctly there. In `cmd.exe`, single quotes do not protect spaces, so the JSON gets split into multiple arguments.

Open PowerShell in your virtual environment and navigate to the example directory:

```powershell
cd ras_to_srh
```

If you must use `cmd.exe`, replace single quotes with outer double quotes and escape the inner double quotes:

```cmd
hmt-cli srh_to_vtk --args "{\"srhhydro_file\": \"Muncie.srhhydro\", \"output_file\": \"my_case_mesh.vtk\", \"mesh_only\": true}"
```

See further instructions in each of the examples.
