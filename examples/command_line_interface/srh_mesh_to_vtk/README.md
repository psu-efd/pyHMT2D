## Readme

This directory contains examples for exporting only the SRH-2D **mesh** to VTK using the `hmt-srh-to-vtk` command with the `--mesh-only` option.

Make sure you are in a Windows terminal with the correct Python virtual environment activated (where *pyHMT2D* is installed). Then, in this example's directory, you can run:

```bash
hmt-srh-to-vtk --mesh-only my_case.srhhydro my_case_mesh.vtk
```

Arguments:

- **srhhydro_file**: SRH-2D case file in `.srhhydro` format (e.g., `my_case.srhhydro`).
- **vtk_file**: Output VTK file name; must include the `.vtk` extension (e.g., `my_case_mesh.vtk`).

For help and version:

```bash
hmt-srh-to-vtk -h
hmt-srh-to-vtk --version
```
