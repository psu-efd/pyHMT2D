## Readme

This directory contains examples for the `hmt-srh-mesh-to-vtk` command of *pyHMT2D*.

Make sure you are in a Windows terminal with the correct Python virtual environment activated (where *pyHMT2D* is installed). Then, in this example's directory, you can run:

```bash
hmt-srh-mesh-to-vtk my_case.srhhydro my_case_mesh.vtk
```

Here, the command `hmt-srh-mesh-to-vtk` takes two positional arguments:

- `srhhydro_file`: SRH-2D case file in .srhhydro format (e.g., `my_case.srhhydro`).
- `vtk_file`: Output VTK file name; must include the .vtk extension (e.g., `my_case_mesh.vtk`).

You can also type the following for more information:

```bash
hmt-srh-mesh-to-vtk -h
hmt-srh-mesh-to-vtk --version
```