---
name: hmt-open
description: Open a pyHMT2D hydraulic model project (SRH-2D or HEC-RAS) and display its materials, boundary conditions, and available result variables. Use this when the user wants to open, load, or inspect a hydraulic model.
---

Open a pyHMT2D hydraulic model project and display its contents.

## Usage
`/hmt-open [project_file]`

The project file path is optional — if omitted, this skill scans the current directory.

## Steps

1. **Find the project file.**
   If the user provided a path in the command, use it directly.
   Otherwise, scan the current directory:
   ```bash
   hmt-cli list_model_files --args '{"directory": "."}'
   ```
   Show the found files as a table (file, model type, case name).
   If multiple files are found, ask the user which one to open.
   If none are found, tell the user no model files were detected and stop.

2. **Open the project.**
   ```bash
   hmt-cli get_project_info --args '{"project_file": "<path>"}'
   ```
   On error: report the message and stop.
   On success: report model type, case name, number of cells, number of material
   zones, and number of boundary conditions.

3. **Show materials.**
   ```bash
   hmt-cli get_materials
   ```
   Present as a Markdown table: ID | Name | Manning's n

4. **Show boundary conditions.** (run after step 3 completes — do NOT run in parallel)
   ```bash
   hmt-cli get_boundary_conditions
   ```
   Present as a Markdown table: ID | Name | Type | Value | Units

5. **Look for a result file.** (run after step 4 completes — do NOT run in parallel)
   Check the same directory for `.h5`, `.hdf`, or `.hdf5` files.
   If found, offer to inspect available variables:
   ```bash
   hmt-cli get_result_variables --args '{"result_file": "<result_file>"}'
   ```
   Report: variable names and number of time steps.

   > **Important**: Each `hmt-cli` call spawns a separate process that opens its own
   > HEC-RAS COM instance. Always run these calls **sequentially** to avoid
   > simultaneous instances that can cause file I/O conflicts.

6. **Suggest next steps** based on what was found:
   - If result file present: "Run `/hmt-results` to query simulation output."
   - If no result file: "Run `/hmt-run` to execute the simulation."
   - To change parameters: "Run `/hmt-modify`."

The session is automatically saved to `.hmt_session.json` by each `hmt-cli` call.
