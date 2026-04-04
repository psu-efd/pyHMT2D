# RAS-2D calibration

This example demonstrates how to calibrate Manning's n for the river channel (zone 2) and
left_2 (zone 4) in the Muncie 2D case using pyHMT2D's `HEC_RAS_Model` class and `RAS_2D_Data`
(accessed via `HEC_RAS_Project`). The calibration uses high water marks from `HWMs.dat`.

The calibration is performed with the `gp_minimize` function from the scikit-optimize library.
The basic steps are:

- Define the objective function: RMSE between simulated WSE and high water marks.
- Define parameter bounds for Manning's n (channel and left_2 zones).
- Define the computational budget (maximum number of objective evaluations).
- Run calibration with `gp_minimize`:
  - The base case is in the directory `base_case`.
  - Copy `base_case` to a temporary run directory for each evaluation.
  - Modify Manning's n in the HEC-RAS geometry HDF using `RAS_2D_Data.modify_ManningsN()`.
  - Run HEC-RAS via `HEC_RAS_Model`.
  - Convert results to VTK using the project-centric API:
    ```python
    from pyHMT2D.Hydraulic_Models_Data.RAS_2D.HEC_RAS_Model import HEC_RAS_Project
    project = HEC_RAS_Project("Muncie2D.prj")
    plan = project.get_plan("p01")
    ras_2d_data = plan.load_results()
    ras_2d_data.saveHEC_RAS2D_results_to_VTK(lastTimeStep=True)
    ```
  - Probe simulated WSE at high water mark locations and compute RMSE.
  - Record parameters and error for each evaluation.
- Save all evaluations to `calibration_results.csv`.
- Copy `base_case` to `calibrated_case` and run with the best parameters found.

**Run the demo**

From this directory (where `base_case` and `HWMs.dat` are located):

```bash
python demo_calibration_RAS_2D.py
```

The script uses pyHMT2D's `HEC_RAS_Model`, `HEC_RAS_Project`, `RAS_2D_Data`, and `vtkHandler`.

## Note
- The calibration is done in serial. It can be parallelized using the `joblib` library together
  with scikit-optimize's `Optimizer` class.
- The calibration script was initially generated with AI assistance (Cursor).
