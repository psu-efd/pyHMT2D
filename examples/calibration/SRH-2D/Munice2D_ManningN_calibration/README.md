# RAS-2D calibration

This example demonstrates how to calibrate the Manning's n for the river channel (zone 2) and left_2 (zone 4) in the Muncie 2D case using pyHMT2D's `RAS_2D_Model` class and `RAS_2D_Data` classes. The calibration is based on high water marks provided in the file `HWMs.dat`.

The calibration is performed using the gp_minimize function from the scikit-optimize library. The parameters to be calibrated are the Manning's n for the river channel (zone 2) and left_2 (zone 4). The basic steps are:    
- Define the objective function to be minimized: error between the high water marks and the predicted high water marks.
- Define the bounds for the parameters to be calibrated.
- Define the computational budget: the maximum number of function evaluations.
- Start the calibration using the gp_minimize function:
  - The initial case is in the directory `base_case`.
  - Modify the Manning's n values in HEC-RAS case files using the `RAS_2D_Model` class.
  - Run the HEC-RAS model to get the results.
  - Convert the HEC-RAS results to VTK format.
  - Probe the simulated water surface elevation at the high water marks and calculate the error between the simulated and the measured high water marks.
  - Return the error as the objective function value.
  - Get the next guess for the parameters to be calibrated from the gp_minimize function. Record the parameters and the objective function value as a function of the number of iterations.
- The results (calibration process) are saved to the file `calibration_results.csv`.
- Make a copy of the base case to a directory named `calibrated_case` and run the case with the calibrated Manning's n values.

**Run the demo**

From this directory (where `base_case` and `HWMs.dat` are located):

```bash
python demo_calibration_RAS_2D.py
```

Requires: `pyHMT2D` and `scikit-optimize` (`pip install scikit-optimize`). The script uses only pyHMT2D’s basic features: `HEC_RAS_Model`, `RAS_2D_Data`, and `vtkHandler` for running HEC-RAS, converting results to VTK, and probing WSE at the high water mark points. 



## Note
The calibration script was generated in Cursor with the following prompt:
```
I need to create a demonstration for the calibration of a HEC-RAS 2D case. The base case has already been prepared in "base_case". The steps and requirements are in @examples/calibration/RAS-2D/Munice2D_ManningN_calibration/README.md The Python script should use pyHMT2D package's basic funcationalities. In the "examples", I already have demonstration cases on how to control HEC-RAS 2D runs, process HEC-RAS 2D results to vtk, sampling at ponts on vtk, etc. 
```
