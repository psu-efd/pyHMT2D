## Readme

This directory contains examples for the `hmt-srh-to-vtk` command of *pyHMT2D*: converting SRH-2D **results** (solution data from an HDF file) to VTK, not mesh-only export.

Make sure you are in a Windows terminal with the correct Python virtual environment activated (where *pyHMT2D* is installed). Then, in this example's directory, you can run:

```bash
hmt-srh-to-vtk my_case.srhhydro my_case_XMDFC.h5
```

Arguments:

- **srhhydro_file**: SRH-2D case file in `.srhhydro` format (e.g., `my_case.srhhydro`).
- **hdf_file**: SRH-2D result file in HDF format (e.g., `my_case_XMDFC.h5`). The command reads this file and exports the solution to VTK.

Optional flags (for result conversion):

- `--nodal`: use nodal (node) data instead of cell-centered (default: False).
- `--all-timesteps`: convert all time steps to VTK; if not set, only the last time step is converted.

To export only the mesh (no result data), use the `--mesh-only` option and pass an output VTK filename as the second argument; see the `srh_mesh_to_vtk` example directory.

For help and version:

```bash
hmt-srh-to-vtk -h
hmt-srh-to-vtk --version
```
