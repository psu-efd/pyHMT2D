---
name: hmt-export
description: Export simulation results and/or the computational mesh to VTK files for ParaView visualization. Use this when the user wants to export, visualize, or create VTK files.
---

Export simulation results and/or the computational mesh to VTK files for ParaView visualization.

## Prerequisites
A project must be open. Results must be loaded (run `/hmt-results` first if needed).

## Steps

1. **Check session state.**
   If `result_loaded` is false in `.hmt_session.json`, load results first:
   ```bash
   hmt-cli read_results --args '{"result_file": "<path>", "timestep": -1}'
   ```

2. **Export results to VTK.**
   ```bash
   hmt-cli export_to_vtk --args '{"timestep": -1, "output_dir": "./vtk_output"}'
   ```
   Report: list of VTK file paths created.

3. **Optionally export the mesh only** (SRH-2D only).
   Ask: "Do you also want to export just the computational mesh (no result data)?"
   If yes:
   ```bash
   hmt-cli export_mesh_to_vtk
   ```

4. **Report all exported file paths.**
   "Open these files in ParaView to visualize the results."
