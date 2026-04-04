---
name: hmt-results
description: Load simulation results and query values at specific points, domain-wide statistics, flood extent, or cross-section profiles. Use this when the user wants to read, query, probe, or analyze simulation results.
---

Load simulation results and query values at points, domain-wide statistics, flood extent, or cross-section profiles.

## Prerequisites
A project must be open. A result file (`.h5` or `.hdf`) must exist.

## Steps

1. **Find the result file.**
   Check `.hmt_session.json` for `result_file`. If not set, inspect any `.h5` / `.hdf` file:
   ```bash
   hmt-cli get_result_variables --args '{"result_file": "<path>"}'
   ```
   If multiple result files exist, ask the user which one to use.

2. **Load the results.**
   ```bash
   hmt-cli read_results --args '{"result_file": "<path>", "timestep": -1}'
   ```
   Report: available variable names, number of time steps.

3. **Ask the user what they want to query.** Offer these options:

   **a) Value at a point** — ask for x, y, variable name:
   ```bash
   hmt-cli get_value_at_point --args '{"x": <x>, "y": <y>, "variable": "<var>"}'
   ```

   **b) Domain-wide statistics** — ask for variable name:
   ```bash
   hmt-cli get_result_statistics --args '{"variable": "<var>"}'
   ```
   Report: min, max, mean, std, number of cells.

   **c) Flood extent** — ask for depth threshold (default 0.01 m):
   ```bash
   hmt-cli get_flood_extent --args '{"depth_threshold": <threshold>}'
   ```
   Report: flooded area in m² and km², percentage of domain flooded.

   **d) Cross-section profile** — ask for two (x,y) endpoints and variable:
   ```bash
   hmt-cli get_cross_section_profile --args '{"x1": <x1>, "y1": <y1>, "x2": <x2>, "y2": <y2>, "variable": "<var>", "n_points": 50}'
   ```
   Present as a table: Distance (m) | X | Y | Value.

4. **Offer further queries** or suggest `/hmt-export` to create VTK files.
