---
name: hmt-modify
description: Modify Manning's n roughness coefficients, inlet flow rates, or exit water surface elevations for the currently open pyHMT2D project and save the changes. Use this when the user wants to change model parameters.
---

Modify Manning's n roughness, inlet flow, or exit water surface elevation for the open project.

## Prerequisites
A project must be open. If not, run `/hmt-open` first.

## Steps

1. **Show current materials and boundary conditions** so the user knows what is available.
   Run sequentially — do NOT run in parallel (each call opens a HEC-RAS instance):
   ```bash
   hmt-cli get_materials
   ```
   ```bash
   hmt-cli get_boundary_conditions
   ```
   Present both as Markdown tables.

2. **Ask the user** what to change. They may specify:
   - Manning's n values by material name or ID and new value
   - Inlet flow by BC ID and new value (SRH-2D only)
   - Exit WSE by BC ID and new value (SRH-2D only)

3. **Apply Manning's n changes** (if any).
   Map material names to IDs from the materials table shown in step 1.
   ```bash
   hmt-cli set_manning_n --args '{"material_ids": [...], "values": [...], "material_names": [...]}'
   ```
   Report the old and new values for each material.

4. **Apply inlet flow changes** (if any, SRH-2D only).
   ```bash
   hmt-cli set_inlet_flow --args '{"bc_ids": [...], "values": [...]}'
   ```
   Note: skip this step for HEC-RAS models.

5. **Apply exit WSE changes** (if any, SRH-2D only).
   ```bash
   hmt-cli set_exit_wse --args '{"bc_ids": [...], "values": [...]}'
   ```

6. **Save the modified input files.**
   ```bash
   hmt-cli save_modified_inputs
   ```
   Report the file path where changes were written.

7. **Confirm** all changes with a summary table: Parameter | Old value | New value.
   Suggest running `/hmt-run` to execute the simulation with the new parameters.
