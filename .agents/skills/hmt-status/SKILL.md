---
name: hmt-status
description: Show the current pyHMT2D session state — which project is open, whether results are loaded, solver configuration, and suggested next steps. Use this when the user asks what is open, what the current state is, or what to do next.
---

Show the current pyHMT2D session state and suggest next steps.

## Steps

1. **Check for an active session.**
   Look for `.hmt_session.json` in the current directory.
   If absent: tell the user there is no active session and suggest `/hmt-open`.

2. **Query the live session.**
   ```bash
   hmt-cli get_simulation_status
   ```

3. **Present a session summary** as a formatted report:
   - Model type (SRH-2D / HEC-RAS)
   - Project file (show basename only for readability)
   - Results loaded (yes / no)
   - Available result variables (comma-separated list)
   - Number of result time steps
   - Last VTK file exported (basename, or "none")
   - Solver configured (yes if `_srh_path` or `_hecras_version` is set, else no)

4. **Suggest next steps** based on current state:
   - Project not open → "Run `/hmt-open` to open a project."
   - Project open, no results → "Run `/hmt-run` to simulate, or `/hmt-results` to load existing results."
   - Results loaded, no VTK → "Run `/hmt-results` to query results or `/hmt-export` to create VTK files."
   - VTK available → "Run `/hmt-results` to probe values, or `/hmt-calibrate` / `/hmt-monte-carlo` for analysis."
   - To close the HEC-RAS GUI without clearing the session: `hmt-cli exit_model`
   - To reset everything and start fresh: `hmt-cli close_session`
