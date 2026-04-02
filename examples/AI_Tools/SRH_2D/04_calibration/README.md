# Example 4: Manning's n Calibration

This example demonstrates how to calibrate the Manning's n for the channel material and the left_2 material (roughness) zones in the SRH-2D model using the AI Tools. For example, you can use the following prompt to get the information you need:

```text
I want to calibrate Manning's n for an SRH-2D model. The config is in ..\hmt_config.json, the original model files are
in base_case\Muncie.srhhydro, and the high water mark observations are in .\HWMs.dat. Please load the config, open the
project, check the observation file format, then calibrate Manning's n for the channel and left_2 material zones. Their respective value ranges are 0.02 to 0.08 for channel and 0.02 to 0.08 for left_2. Use a computing budget of 11 iterations. Save the calibration history to calibration_history.csv.
```

