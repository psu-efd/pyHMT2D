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
Efficient and accurate flood simulation is increasingly important to human society as climate keeps changing and floods of all magnitudes seem to occur more frequent. Two-dimensional (2D) hydraulic models, replacing one-dimensional (1D) models, have become the work horse for most academic research and engineering purposes in practice. Many agencies, such as U.S. DOT, Bureau of Reclamation (USBR), FEMA, and U.S. Army Corp of Engineers (USACE), have developed and promoted 2D hydraulic models to fulfill their respective missions. Example 2D models are SRH-2D (USBR) and HEC-RAS 2D (USACE). Over the past several decades, these models have been greatly improved by their respective agencies independently. However, limitations still exist in these individual models. For example, it is desirable to have fair and meaningful model result inter-comparison. But due to the differences in numerical schemes and model approximations, such comparison is currently difficult, if not impossible. Another example of limitation is that most of these models evolved from legacy codes. They lack the flexibility to cope with the rapidly evolving landscape of computational science and utilize the power of modern computational ecosystems such as Python and R. 

`pyHMT2D` is a Python package developed to partially remove these limitations. It stands for "2D Hydraulic Modeling Tools with Python". Objected-oriented programming is utilized to encapsulate the hydraulic models and simulation data (Fig. 1). A user only needs to interface with `pyHMT2D` using the high-level APIs implemented in a handful of classes in the package. In addition to these core classes, `pyHMT2D` also provides a suite of miscellaneous tools to perform various tasks essential to hydraulic modeling tasks, such as creation and manipulation of geo-referenced terrain data. To compare model results and provide an alternative visualization means, VTK [@vtkbook2006] is used as a common format. `pyHMT2D` can read results from all supported hydraulic models and translates them to VTK. Results in VTK format can be opened in ParaView [@Ayachit2015] for visualization and further analysis. Another feature that makes `pyHMT2D`more attractive is the use of Python as a "glue" language to expose all supported hydraulic models to the large amount of powerful optimization packages. It can facilitate automatic and more intelligent model calibrations. In summary, `pyHMT2D` can be used to control and (semi)automate 2D hydraulic modeling, pre-/post-processing simulation results, convert simulation results into the unified format of VTK [@vtkbook2006] and automatic calibration.

![Code scheme diagram for the `pyHMT2D` package.](./figures/package_scheme_diagram.png)


# Statement of need

Floods are one of the most common natural hazards which affect hundreds of millions of people worldwide and cost billions of dollars in damage. In the United States, about 40\% of recorded 35,000 disasters since 1900 were major floods and related storms [@Cigler2017]. With
climate change, floods are becoming more frequent and severe. To better understand the dynamics of floods and mitigate their impact, hydraulic models have been used extensively. In engineering practice, hydraulic models are also widely used to design and evaluate infrastructure, such as bridges and roads, against the damaging effect of flowing water. Flow in rivers and streams is naturally three-dimensional. However, it is computationally prohibitive to simulate three-dimensional flow field for large scale applications.   Two-dimensional (2D) hydraulic modeling is the prevailing approach for both engineering practice and research. Many 2D models exist and the fundamental physical governing equations implemented are very similar. Indeed, the depth-averaged shallow water equations are the backbone of most 2D hydraulic simulation models [@Liu2008]. However, different 2D hydraulic models have their own unique numerical scheme, mesh and result format, and even some variations to the governing equations. The diverse, isolated, and most time incompatible ecosystems surrounding these models can create problems. 

Due to the importance of 2D hydraulic models to the society for flood prediction, prevention, mitigation, and engineering design of infrastructure, there is a great need to make them comparable. In addition, as for all other computational models, there is a persistent need to make them more efficient, easy to use, automatic to be calibrated. Therefore, the motivations of developing `pyHMT2D` are as follows:

- One major motivation is to efficiently and automatically run 2D hydraulic modeling simulations, for example, batch simulations to intelligently and efficiently calibrate models. Many of the 2D models have some automation to certain degree. However, these models and their GUIs are closed source. Therefore, a modeler is limited to what he/she can do. Automated calibration using the Python code in `pyHMT2D` also makes it possible for these models to better ingest vast amount of data generated by millions of sensors deployed in the field.
- `pyHMT2D` also serves as a bridge between 2D hydraulic models and the Python universe where many powerful libraries exist, for example statistics, optimization, machine learning, GIS, and parallel computing.
- `pyHMT2D` can make model inter-comparison and evaluation an easy task. Almost all 2D hydraulic models solve the shallow-water equations. However, every model does it differently. How these differences manifest in their results and how to quantify/interpret the differences are of great interest to practitioners.
- The read/write and transformation of 2D hydraulic model results can be used to feed other models which use the simulated flow field, for example external water quality models and fish migration models.
- Most 2D models have good user interface and they have the capability to produce good result visualizations and analysis. However, with `pyHMT2D` and the power of the VTK library, 2D hydraulic modeling results can be visualized and analyzed with more flexibility and efficiency.

Previously, tools have been developed to control hydraulic models and process their simulation results. Most of these existing tools use similar scripting languages, such as VBA, Matlab, and Python. For example, @MoyaQuirogaEtAl2013 showed the evaluation of HEC-RAS modeling uncertainties with the combination of Delphi, VB.NET, and Windows command line batch scripts. In @Goodell2014 and its companion code, VBA was used to interact with HEC-RAS and a variety of demonstration examples were provided. @LeonGoodell2016 developed several MATLAB scripts for controlling HEC-RAS, which include input file reading and writing, output file reading, plotting, and parallel computation. Similar implementation using Python has been reported in @Dysarz2018. If one searches on GitHub with the keyword "HEC-RAS", there are many repositories and code snippets which also implement certain functionality mentioned above and they are for one particular model only. However, `pyHMT2D` is designed to be more general and extensible. The class hierarchy shown in Fig. 1 can easily add other hydraulic models into `pyHMT2D`. Model inter-comparison functionality is also missing in existing tools.  

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

