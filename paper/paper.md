---
title: 'pyHMT2D: Python-based two-dimensional hydraulic modeling tools'
tags:
  - Python
  - rivers
  - flood
  - SRH-2D
  - HEC-RAS
  - computational hydraulics
authors:
  - name: Xiaofeng Liu^[Corresponding author]
    orcid: 0000-0002-8296-7076
    affiliation: "1"
affiliations:
 - name: Department of Civil and Environmental Engineering, Institute of Computational and Data Sciences, The Pennsylvania State University
   index: 1
date: 11 April 2021
bibliography: paper.bib
---

# Summary
Efficient and accurate flood simulation is increasingly important to human society as climate keeps changing and floods of all magnitudes seem to occur more frequent. Two-dimensional (2D) hydraulic models, replacing one-dimensional (1D) models, have become the work horse for most academic research and engineering purposes in practice. Many agencies, such as U.S. DOT, Bureau of Reclamation (USBR), FEMA, and U.S. Army Corp of Engineers (USACE), have developed and promoted 2D hydraulic models to fulfill their respective missions. Example 2D models are SRH-2D (USBR) and HEC-RAS 2D (USACE). Over the past several decades, these models have been greatly improved by their respective agencies independently. However, limitations still exist in these individual models. For example, it is desirable to have fair and meaningful model result inter-comparison. But due to the differences in numerical schemes and model approximations, it is currently impossible. Another example of limitation is that most of these models evolved from legacy codes. They lack the flexibility to cope with the rapidly evolving landscape of computational science and utilize the power of modern computational ecosystems such as Python and R. 

`pyHMT2D` is a Python package developed to partially remove these limitations. It stands for "2D Hydraulic Modeling Tools with Python". Objected-oriented programming is utilized to encapsulate the hydraulic models and simulation data. A user only needs to interface with `pyHMT2D` using the high-level APIs implemented in a handful of classes in the package. In addition to these core classes, `pyHMT2D` also provides a suite of miscellaneous tools to perform various tasks essential to hydraulic modeling tasks, such as creation and manipulation of geo-referenced terrain data. To compare model results and provide an alternative visualization tool, VTK [@vtkbook2006] is used as a common format. `pyHMT2D` can read results from all supported hydraulic models and translates them to VTK. Another feature that makes `pyHMT2D`more attractive is the use of Python as a "glue" language to expose all supported hydraulic models to the large amount of powerful optimization packages. It can facilitate automatic and more intelligent model calibrations. In summary, `pyHMT2D` can be used to control and (semi)automate 2D hydraulic modeling, pre-/post-processing simulation results, convert simulation results into the unified format of VTK [@vtkbook2006] and automatic calibration.


# Statement of need


SWEs are the backbone of most flood simulation models [@Liu2008].

# Functionality

`pyHMT2D` requires the user to have some basic knowledge about Python programming, which is 


### Basic Functionality


### Advanced Functionality


# Dependencies

`pyHMT2D`'s  core functionality relies the following Python packages: h5py [@collette2013], VTK [@vtkbook2006], GDAL [@gdal2020], NumPy [@harris2020], and Matplotlib [@hunter2007]. Optionally, pywin32 [@hammond2000] is needed to control the simulation of HEC-RAS and affine [@Gillies2014] is needed to read HEC-RAS results. 

# Acknowledgements

I thank Dr. Yong Lai at U.S. Bureaul of Reclamation and Scott Hogan at U.S. Federal Highway Administration for their valuable knowledge and suggestions regarding SRH-2D. I also thank Gary Brunner from U.S. Army Corps of Engineers for his time answering my questions regarding HEC-RAS.

`pyHMT2D` was developed as a tool box for several past and ongoing research projects in our group. Though `pyHMT2D` did not directly receive financial support from these projects nor it is in the workscope, `pyHMT2D` did benefit from the knowledge on various aspects of 2D hydraulic models accumulated through these projects. 

# References

