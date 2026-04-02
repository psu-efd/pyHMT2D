# Example 3: Query Results

This example demonstrates how to query the results of an SRH-2D project using the AI Tools. For example, you can use the following prompt to get the information you need:

```text
The config file is ..\hmt_config.json and the SRH-2D model is Muncie.srhhydro in the current directory. Load the config, open the project, read the results from the XMDFC file Muncie_XMDFC.h5, print the list of available result variables, calculate the domain-wide statistics for the water surface elevation and depth, get the water surface elevation at a specific point (e.g., a bridge pier location at x=411100.0, y=1803450.0), calculate the flood extent (area above a depth threshold) above a depth threshold (e.g., 0.1 ft), and extract a cross-section profile (from x=410500.0, y=1803000.0 to x=411500.0, y=1803000.0). Finally, export the results to VTK for ParaView visualisation.
```