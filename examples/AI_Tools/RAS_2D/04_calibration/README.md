# Example 4: Manning's n Calibration

This example demonstrates how to calibrate the Manning's n for the channel material and the left_2 material (roughness) zones in a HEC-RAS 2D model using the AI Tools. For example, you can use the following prompt to get the information you need:

```text
I want to calibrate Manning's n for a HEC-RAS 2D model. The original model files are
in "base_case" directory , and the high water mark observations are in "HWMs.dat". Please check the observation file format, then calibrate Manning's n for the channel and left_2 material zones. Their respective value ranges are 0.02 to 0.08 for channel and 0.02 to 0.08 for left_2. Use a computing budget of 11 iterations. Save the calibration history to calibration_history.csv.
```

