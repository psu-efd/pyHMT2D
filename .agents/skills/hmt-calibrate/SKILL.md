---
name: hmt-calibrate
description: Run automated Manning's n calibration against observed water surface elevations using Bayesian optimization or Nelder-Mead. Use this when the user wants to calibrate, optimize, or fit model parameters to observations.
---

Run automated Manning's n calibration against observed water surface elevations.

## Prerequisites
- Project opened from a `base_case/` subdirectory
- Observation CSV: columns `x, y, wse` or `name, x, y, wse` (comment lines start with `#`)
- Solver paths configured in `hmt_config.json`

## Steps

1. **Check that a project is open.**
   If not, ask the user to run `/hmt-open` pointing to `base_case/<project_file>`.

2. **Validate the observation file.**
   Ask the user for the observation CSV path (e.g., `HWMs.dat`).
   ```bash
   hmt-cli check_observation_format --args '{"csv_file": "<path>"}'
   ```
   Show the preview table. On error: report and ask for a corrected path.

3. **Show available materials.**
   ```bash
   hmt-cli get_materials
   ```

4. **Ask the user for calibration parameters.** For each material to calibrate:
   - Material name, minimum n, maximum n, optional initial guess

5. **Build validated parameter specifications.**
   ```bash
   hmt-cli build_param_specs --args '{"specs": [
     {"type": "manning_n", "material_name": "<name>", "min": <min>, "max": <max>, "initial": <init>},
     ...
   ]}'
   ```
   Show the specs and confirm with the user.

6. **Optional: single test evaluation** at initial values to verify the setup.
   ```bash
   hmt-cli evaluate_parameters --args '{"param_specs": <specs_with_initial_values>, "observation_csv": "<path>"}'
   ```
   If RMSE = 1e6, the simulation failed — check `pyHMT2D.log` before proceeding.

7. **Ask for settings:** number of iterations (default: automatic) and method (`"gp"` or `"nelder-mead"`).

8. **Run the calibration.**
   ```bash
   hmt-cli run_calibration --args '{"param_specs": <specs>, "observation_csv": "<path>", "n_iterations": <N>, "method": "gp"}'
   ```
   Periodically show progress while waiting:
   ```bash
   tail -20 calib_progress.log
   ```

9. **Report results:** best parameter values table, best RMSE, iterations completed, history CSV path.

**Troubleshooting:**
- RMSE = 1,000,000 → check `pyHMT2D.log` for solver errors
- Material not found → name must match exactly; copy from `hmt-cli get_materials` output
