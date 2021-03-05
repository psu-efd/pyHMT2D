
[comment]: <> (
<p align="center">
  <a href="https://github.com/psu-efd/pyHMT2D"><img alt="pyHMT2D" src="logo/pyhmt2d_logo_vel_color_vector.png" width="60%"></a>
</p>
)

![image](logo/pyhmt2d_logo_vel_color_vector.png)
<p align="center">2D Hydraulic Modeling Tools in Python.</p>

# pyHMT2D - <ins>Py</ins>thon <ins>H</ins>ydraulic <ins>M</ins>odeling <ins>T</ins>ools - 2D

pyHMT2D is a Python package developed to control and (semi)automate 2D hydraulic modeling, and pre-/postprocessing simulation results. Currently, the following 2D hydraulic models are supported:
- SRH-2D
- HEC-RAS 2D

## Motivations

Two-dimensional (2D) hydraulic modeling, replacing one-dimensional (1D) modeling, has become the work horse for most engineering purposes in practice. Many agencies, such as U.S. DOT, Bureau of Reclamation (USBR), FEMA, and U.S. Army Corp of Engineers (USACE), have developed and promoted 2D hydraulic models to fulfill their respective missions. Example 2D models are SRH-2D (USBR) and HEC-RAS 2D (USACE). The motivations of this package are as follows:

- One major motivation of this package is to efficiently and automatically run 2D hydraulic modeling simulations, for example, batch simulations to calibrate model runs. Many of the 2D models have some automation to certain degree. However, these models and their GUIs are closed source. Therefore, a modeler is limited to what he/she can do. 
- Most 2D models have good user interface and they have capability to produce good result visualizations and analysis. However, with this package and the power of the VTK library, 2D hydraulic modeling results can be visualized and analyzed with more flexibilty and efficiency. 
- This package also serves as a bridge between 2D hydraulic models and the Python universe where many powerful libaries exist, for example statistics, machine learning, GIS, and parallel computing.
- The read/write and tranformation of 2D hydraulic model results can be used to feed other models which use the simulated flow field, for example external water quality models and fish models.
- Model inter-comparison and evaluation. Almost all 2D hydraulic models solve the shallow-water equations. However, every model does it differently. How these differences manifest in their results and how to quantify/interpret the differences are of great interest to practitioners.

## Features

For SRH-2D modeling:
- read SRH-2D results if they are not in VTK format
- convert SRH-2D results to VTK format
- sample and probe simulation results (with the functionalities of VTK library)
- control SRH-2D simulations (TODO)

For HEC-RAS 2D modeling:
- read RAS2D results (HDF files)
- convert RAS2D results to VTK, one of the most popular format for scientific data
  - point and cell center data (depht, water surface elevation, velocity, etc.)
  - interpolate between point and cell data
  - face data (e.g., subterrain data)
- convert RAS2D mesh, boundary conditions, and Manning's n data into SRH-2D format such that SRH-2D and HEC-RAS 2D can run a case with exactly the same mesh for comparison purpose. Additionally, HEC-RAS 2D can be used as a mesh generator for SRH-2D. 
- sample and probe simulation results (with the functionalities of VTK library)
- control HEC-RAS 2D simulations (TODO)

With the control and automation capability above, it is much easier to do the following:
- automatic calibration of models with any available optimzation and calibration Python packages
- Monte-Carlo simulations with scripting and Python's statistic libraries
- ...

Other features:
- calculate the difference between simulation results (if they are on the same mesh)

## Requirements

This package uses the following libraries:
- [meshio](https://github.com/nschloe/meshio)
- [h5py](https://www.h5py.org/)
- [vtk](https://github.com/Kitware/VTK)

## Installation

## Example Usage

More examples can be found in the "test" directory.

## Limitations

For SRH-2D:
- This package is developed and tested with SRH-2D v3.3; other versions may work but has not been tested.

For HEC-RAS 2D:
- Only one 2D flow area is supported.
- Only 2D flow area information is processed; others such as 1D channels and structures are ignored.
- Only flow data is processes; others such as sediment and water quality are ignored.
- This package is developed and tested with HEC-RAS v5.0.7 and v6.0 beta; other versions may work but has not been tested.

## License

MIT


## Author

Xiaofeng Liu, Ph.D., P.E.  
Associate Professor

Department of Civil and Environmental Engineering  
Institute of Computational and Data Sciences  
Penn State University

223B Sackett Building, University Park, PA 16802

Web: http://water.engr.psu.edu/liu/

## Contributors
(To be added)
