---
name: hmt-convert
description: Convert between hydraulic model formats. Supports HEC-RAS 2D to SRH-2D mesh conversion, and SRH-2D results or mesh to VTK. These tools are stateless — no open project is required. Use this when the user wants to convert formats, e.g. "convert RAS to SRH", "export SRH-2D to VTK", or "convert mesh to VTK".
---

Convert between hydraulic model formats using stateless conversion tools.
No open project or session is required.

## Usage
`/hmt-convert [description of what to convert]`

## Supported conversions

| Conversion | Tool |
|---|---|
| HEC-RAS 2D mesh + Manning's n → SRH-2D (`.srhgeom`, `.srhmat`) | `ras_to_srh` |
| SRH-2D HDF results → VTK | `srh_to_vtk` |
| SRH-2D mesh only → VTK | `srh_to_vtk` with `mesh_only: true` |

## Steps

1. **Determine the conversion type** from the user's request or arguments.
   If unclear, ask:
   - "Do you want to convert HEC-RAS to SRH-2D, or export SRH-2D results/mesh to VTK?"

2. **Collect required file paths.**
   If the user did not provide them, scan the current directory for relevant files:
   ```bash
   hmt-cli list_model_files --args '{"directory": "."}'
   ```
   For `ras_to_srh`: need `.hdf` result file, terrain `.tif`, and desired output case name.
   For `srh_to_vtk`: need `.srhhydro` file and either the HDF result file or output VTK filename.

3. **Run the conversion.**

   **HEC-RAS → SRH-2D:**
   ```bash
   hmt-cli ras_to_srh --args '{"ras_hdf_file": "<path>", "terrain_tif_file": "<path>", "srh_case_name": "<name>"}'
   ```
   Report the generated `.srhgeom` and `.srhmat` files.
   Remind the user that a `.srhhydro` boundary-condition file must be created separately.

   **SRH-2D results → VTK (last time step):**
   ```bash
   hmt-cli srh_to_vtk --args '{"srhhydro_file": "<path>", "output_file": "<result.h5>"}'
   ```

   **SRH-2D results → VTK (all time steps):**
   ```bash
   hmt-cli srh_to_vtk --args '{"srhhydro_file": "<path>", "output_file": "<result.h5>", "all_timesteps": true}'
   ```

   **SRH-2D mesh only → VTK:**
   ```bash
   hmt-cli srh_to_vtk --args '{"srhhydro_file": "<path>", "output_file": "<mesh.vtk>", "mesh_only": true}'
   ```

4. **Report output files** and suggest next steps:
   - For HEC-RAS → SRH-2D: "Create a `.srhhydro` file to complete the SRH-2D case, then run `/hmt-open` to open it."
   - For VTK output: "Open the VTK file(s) in ParaView to visualize the results."
