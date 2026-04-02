# Example 5: Monte Carlo Simulation

This example demonstrates how to perform a Monte Carlo simulation for the SRH-2D model using the AI Tools. For example, you can use the following prompt to get the information you need:

```text
The config file is ..\hmt_config.json and the SRH-2D model is .\base_case\Muncie.srhhydro in the current directory. Load the config, open the project, build the parameter specifications for the channel material, generate 100 samples from a truncated-normal distribution with a mean of 0.04 and a standard deviation of 0.005, and a minimum of 0.03 and a maximum of 0.05, run the Monte Carlo simulation, and report the results. For the result analysis, I have a bridge pier monitoring point at x=411100.0, y=1803450.0, and I want to know the probability of the water surface elevation at the exceedance probability of 99%, 90%, 50%, 10%, and 1%. For the whole domain, I want to know the probability of the water depth for each computational cell in the mesh at the exceedance probability of 99%, 90%, 50%, 10%, and 1%; save the results to a VTK file.
```

If you already finished the MC simulations and only want to do the result analysis, you can use the following prompt to get the information you need:

```text
The config file is ..\hmt_config.json. This folder is for the monte carlo simulations. The base case is in base_case\Muncie.srhhydro. I already did the MC simulations and the results are in the vtk files in mc_runs. Now I want to do some postprocessing to get MC statistics. specifically, for a bridge pier monitoring point at x=411100.0, y=1803450.0,  I want to know the probability of the water surface elevation at the exceedance probability of 99%, 90%, 50%, 10%, and 1%. For the whole domain, I want to know the probability of the water depth for each computational cell in the mesh at the exceedance probability of 99%, 90%, 50%, 10%, and 1%; save the results to a VTK file.
```