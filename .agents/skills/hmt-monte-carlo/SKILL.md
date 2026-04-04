---
name: hmt-monte-carlo
description: Run Monte Carlo uncertainty analysis on a hydraulic model by sampling parameters from statistical distributions and computing exceedance probabilities. Use this when the user wants uncertainty analysis, Monte Carlo simulation, or probabilistic results.
---

Run Monte Carlo uncertainty analysis on a hydraulic model.

## Prerequisites
- Project in a `base_case/` subdirectory
- Parameter uncertainty distributions defined by user
- Solver paths configured in `hmt_config.json`

## Steps

1. **Check that a project is open.**
   If not, ask the user to run `/hmt-open` pointing to `base_case/<project_file>`.

2. **Show available materials.**
   ```bash
   hmt-cli get_materials
   ```

3. **Ask the user for uncertain parameters.** For each:
   - Material name or BC ID
   - Distribution: `"truncated_normal"` (default) or `"uniform"`
   - For `truncated_normal`: mean, std, min, max
   - For `uniform`: min, max

4. **Build validated parameter specifications.**
   ```bash
   hmt-cli build_param_specs --args '{"specs": [
     {"type": "manning_n", "material_name": "<name>", "distribution": "truncated_normal",
      "mean": <mean>, "std": <std>, "min": <min>, "max": <max>},
     ...
   ]}'
   ```
   Confirm the distributions with the user.

5. **Ask for MC settings:** number of samples (50–200), random seed (default 42),
   output directory (default `./mc_runs`), parallel processes (default 1).

6. **Generate samples.**
   ```bash
   hmt-cli generate_mc_samples --args '{"param_specs": <specs>, "n_samples": <N>, "random_seed": <seed>, "output_csv": "mc_samples.csv"}'
   ```
   Show a preview of the first 5 sample rows.

7. **Run Monte Carlo simulations.**
   ```bash
   hmt-cli run_monte_carlo --args '{"base_case_dir": "./base_case", "param_specs": <specs>, "n_samples": <N>, "n_processes": <procs>, "random_seed": <seed>, "sample_csv": "mc_samples.csv", "delete_cases": true, "output_dir": "<output_dir>"}'
   ```
   Periodically show progress while waiting:
   ```bash
   tail -20 mc_progress.log
   ```

8. **Report:** successful/failed runs, results JSON path.

9. **Compute statistics.** Ask for observation point coordinates if the user wants point statistics.
   ```bash
   hmt-cli get_mc_statistics --args '{"results_json": "<path>", "observation_points": [{"name": "<name>", "x": <x>, "y": <y>}], "exceedance_probabilities": [99, 90, 50, 10, 1]}'
   ```
   Present exceedance table: Point | P99 | P90 | P50 | P10 | P1

10. **Report** the spatial exceedance VTK path for ParaView visualization.

**Troubleshooting:**
- Many failed cases → check `pyHMT2D.log`
- Variable detection fails → pass `"variable": "<name>"` explicitly
