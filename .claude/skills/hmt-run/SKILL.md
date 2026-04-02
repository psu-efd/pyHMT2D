---
name: hmt-run
description: Run the SRH-2D preprocessor and hydraulic simulation for the currently open project. Use this when the user wants to run, execute, or simulate the model.
---

Run preprocessing (SRH-2D) and the hydraulic simulation for the open project.

## Prerequisites
A project must be open (`/hmt-open`) and input files must be saved.

## Steps

1. **Check the session.** Read `.hmt_session.json` for `model_type`.
   If no project is open, ask the user to run `/hmt-open` first.

2. **Run preprocessing (SRH-2D only).**
   Skip this step for HEC-RAS models.
   ```bash
   hmt-cli run_preprocessing
   ```
   Report elapsed time on success. On error: show the message and stop.

3. **Run the simulation.**
   ```bash
   hmt-cli run_simulation --args '{"timeout_minutes": 120}'
   ```
   For HEC-RAS, pass `"faceless": false` by default. If the user does not want to see the GUI window, pass `"faceless": true`.
   Warn the user this may take several minutes. Report elapsed time and status.
   On error: show the message and suggest checking `pyHMT2D.log`.

4. **Confirm result files were created.**
   ```bash
   hmt-cli list_model_files --args '{"directory": "."}'
   ```

5. **Suggest next steps:**
   "Run `/hmt-results` to query output, or `/hmt-export` to create VTK files for ParaView."
   For HEC-RAS runs with `faceless=false`: "Run `hmt-cli exit_model` to close the HEC-RAS GUI, or `hmt-cli close_session` to reset the session entirely."
