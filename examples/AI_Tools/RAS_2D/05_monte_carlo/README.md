# Example 5: Monte Carlo Simulation

This example demonstrates how to perform a Monte Carlo simulation for the HEC-RAS 2D model using the AI Tools. For example, you can use the following prompt to get the information you need:

```text
I want to do a Monte Carlo simulation for a HEC-RAS case. The base case is in "base_case". Build the parameter specifications for the channel material zone, generate 100 samples from a truncated-normal distribution with a mean of 0.04 and a standard deviation of 0.005, and a minimum of 0.03 and a maximum of 0.05, run the Monte Carlo simulation, and report the results. For the result analysis, I have a bridge pier monitoring point at x=411100.0, y=1803450.0, and I want to know the probability of the water surface elevation at the exceedance probability of 99%, 90%, 50%, 10%, and 1%. For the whole domain, I want to know the probability of the water depth for each computational cell in the mesh at the exceedance probability of 99%, 90%, 50%, 10%, and 1%; save the results to a VTK file.
```

Or if you want to give a more specific, step-by-step instruction, you can use the following:
```text
I want to run a Monte Carlo uncertainty analysis on a HEC-RAS 2D model. The original
  model files are in the base_case directory. Please do the following:                                                                                         1. Load the configuration from ../hmt_config.json.                                      
  2. The base project is in base_case with Muncie2D.prj as the project file.
  3. Set up one uncertain parameter: Manning's n for the channel material zone, using a   
  truncated normal distribution with mean=0.04, std=0.005, min=0.03, max=0.05.
  4. Generate 100 Monte Carlo samples with random seed 42. Save the sample CSV to
  mc_runs/mc_samples.csv.
  5. Run the Monte Carlo ensemble using those samples: base case in base_case/, output to 
  mc_runs/, 1 process (serial), delete individual case directories after each run to save 
  disk space.
  6. Compute exceedance probability statistics on the results using the variable
  Water_Depth_ft, at probe point "probe_A" located at (x=411100.0, y=1803450.0), for      
  exceedance probabilities 99%, 90%, 50%, 10%, and 1%. Save the spatial exceedance VTK to 
  mc_runs/mc_exceedance.vtk.
  7. Exit the HEC-RAS model to release the COM instance.
```

If you already finished the MC simulations and only want to do the result analysis, you can use the following prompt to get the information you need:

```text
This folder is for the monte carlo simulations. The base case is in the base_case directory. I already did the MC simulations and the results are in the vtk files in mc_runs. Now I want to do some postprocessing to get MC statistics. specifically, for a bridge pier monitoring point at x=411100.0, y=1803450.0,  I want to know the probability of the water surface elevation at the exceedance probability of 99%, 90%, 50%, 10%, and 1%. For the whole domain, I want to know the probability of the water depth for each computational cell in the mesh at the exceedance probability of 99%, 90%, 50%, 10%, and 1%; save the results to a VTK file.
```

