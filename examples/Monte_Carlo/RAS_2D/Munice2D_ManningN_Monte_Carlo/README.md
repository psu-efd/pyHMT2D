## Monte Carlo simulations with `HEC_RAS_Model`

This example demonstrates how to use *pyHMT2D*'s `HEC_RAS_Model` and `HEC_RAS_Project` classes
to perform Monte Carlo simulations with HEC-RAS.

In this case, Manning's n for the river channel (zone 2) is sampled from a truncated normal
distribution and varied across simulations.

Each Monte Carlo case:
1. Copies the `base_case` directory to an isolated case directory.
2. Modifies Manning's n in the HEC-RAS geometry HDF using `RAS_2D_Data.modify_ManningsN()`.
3. Runs HEC-RAS via `HEC_RAS_Model`.
4. Converts results to VTK using the project-centric API:
   ```python
   from pyHMT2D.Hydraulic_Models_Data.RAS_2D.HEC_RAS_Model import HEC_RAS_Project
   project = HEC_RAS_Project("Muncie2D.prj")
   plan = project.get_plan("p01")
   ras_2d_data = plan.load_results()
   ras_2d_data.saveHEC_RAS2D_results_to_VTK(lastTimeStep=True)
   ```
5. Copies the VTK result to the `cases/` directory and (optionally) deletes the case directory.

After all simulations are done, the results are post-processed to:
- Calculate WSE exceedance probability at a monitoring point.
- Generate water depth exceedance of probability for all mesh cells and save to
  `water_depth_exceedance_of_probabilities.vtk` (contains 99%, 90%, 50%, 10%, and 1% intervals).

**Run the demo**

Prepare a parameter file with one Manning's n value per line (comma-separated), then:

```bash
python demo_HEC_RAS_Monte_Carlo.py <parameter_file>
```

The simulation can run in parallel (set `bParallel = True` in the script) using
`multiprocessing.Pool`. For parallel runs, configure HEC-RAS to use a single core per case
(*Unsteady Flow Analysis → Options → Computation Options and Tolerances → 1 core*).
