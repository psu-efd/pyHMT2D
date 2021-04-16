.. image:: https://github.com/psu-efd/pyHMT2D/raw/main/docs/source/_static/images/pyhmt2d_logo_vel_color_vector.png
   :target: https://github.com/psu-efd/pyHMT2D
   :alt: pyHMT2D

*pyHMT2D* - Python-based Hydraulic Modeling Tools - 2D
=======================================================

*pyHMT2D* is a Python package developed to control and (semi)automate 2D
hydraulic modeling, and pre-/post-process simulation results.
Currently, the following hydraulic models are supported:

-  `SRH-2D <https://www.usbr.gov/tsc/techreferences/computer%20software/models/srh2d/index.html>`__
-  `HEC-RAS <https://www.hec.usace.army.mil/software/hec-ras/>`__
-  Backwater-1D (a simple toy model for demonstration purpose)

In the future, support for more 2D models may be added.

Motivations
-----------

Two-dimensional (2D) hydraulic modeling, replacing one-dimensional (1D)
modeling, has become the work horse for most engineering purposes in
practice. Many agencies, such as U.S. DOT, Bureau of Reclamation (USBR),
FEMA, and U.S. Army Corp of Engineers (USACE), have developed and
promoted 2D hydraulic models to fulfill their respective missions.
Example 2D models are SRH-2D (USBR) and HEC-RAS 2D (USACE). The
motivations of this package are as follows:

-  One major motivation of this package is to efficiently and
   automatically run 2D hydraulic modeling simulations, for example,
   batch simulations to calibrate model runs. Many of the 2D models have
   some automation to certain degree. However, these models and their
   GUIs are closed source. Therefore, a modeler is limited to what
   he/she can do.
-  Most 2D models have good user interface and they have capability to
   produce good result visualizations and analysis. However, with this
   package and the power of the VTK library, 2D hydraulic modeling
   results can be visualized and analyzed with more flexibility and
   efficiency.
-  This package also serves as a bridge between 2D hydraulic models and
   the Python universe where many powerful libraries exist, for example
   statistics, machine learning, GIS, and parallel computing.
-  The read/write and transformation of 2D hydraulic model results can be
   used to feed other models which use the simulated flow field, for
   example external water quality models and fish models.
-  Model inter-comparison and evaluation. Almost all 2D hydraulic models
   solve the shallow-water equations. However, every model does it
   differently. How these differences manifest in their results and how
   to quantify/interpret the differences are of great interest to
   practitioners.

Features
--------

For SRH-2D modeling:

-  read SRH-2D results (mainly into Python's Numpy arrays)
-  convert SRH-2D results to VTK format, one of the most popular formats for
   scientific data
-  sample and probe simulation results (with the functionalities of VTK
   library)
-  control and modify SRH-2D simulations

For HEC-RAS 2D modeling:

-  read RAS2D results (HDF files, mainly into Python's Numpy arrays)
-  convert RAS2D results to VTK

   -  point and cell center data (depht, water surface elevation,
      velocity, etc.)
   -  interpolate between point and cell data
   -  face data (e.g., subterrain data)

-  convert RAS2D mesh, boundary conditions, and Manning’s n data into
   SRH-2D format such that SRH-2D and HEC-RAS 2D can run a case with
   exactly the same mesh for comparison purpose. As a side bonus, HEC-RAS
   2D can be used as a mesh generator for SRH-2D.
-  sample and probe simulation results (with the functionalities of VTK
   library)
-  control and modify HEC-RAS 2D simulations

With the control and automation capability above, it is much easier to
do the following:

-  automatic calibration of models with available optimization and
   calibration Python packages. Currently, *scipy*'s *optimize* module is supported, which
   includes many local and global optimization methods.
-  Monte-Carlo simulations with scripting and Python’s statistic
   libraries
-  ...

Other features:

-  calculate the difference between simulation results (regardless they are on the same mesh or not)
-  create and manipulate georeferenced terrain data for 2D modeling
-  conversion of 2D model mesh and result to 3D through extrusion (one layer or multiple layers) and
   VTK interpolation. This feature is useful to use 2D simulation result in 3D applications, e.g., fish
   passage design or use 2D result as initial condition for 3D CFD simulations. Currently, conversion
   to `OpenFOAM <https://www.openfoam.com/>`__ is supported through `Gmsh <https://gmsh.info/>`__'s MSH file format.

Requirements
------------

This package uses the following libraries:

-  `h5py <https://www.h5py.org/>`__
-  `vtk <https://github.com/Kitware/VTK>`__
-  `gdal <https://pypi.org/project/GDAL/>`__
-  `pywin32 <https://pypi.org/project/pywin32/>`__ (optional; only if you want to use *pyHMT2D* to control HEC-RAS)
-  `affine <https://pypi.org/project/affine/>`__ (optional; only if you want to use *pyHMT2D* to read HEC-RAS 2D results)

See *pyHMT2D*'s User Manual for how to install these libraries.

Installation
------------

There are several ways to install *pyHMT2D*.

- Install from `pip`

.. code-block:: bash

   $ pip install pyHMT2D

- Directly install from GitHub with `pip`

.. code-block:: bash

   $ pip install git+https://github.com/psu-efd/pyHMT2D.git

- Clone the GitHub repository to your local machine and
  add the local *pyHMT2D*'s directory to your Python path

.. code-block:: bash

   $ git clone https://github.com/psu-efd/pyHMT2D.git

..

     If you use this approach to install *pyHMT2D*, in order to make the Python Interpreter aware of *pyHMT2D*, you need to add the path to *pyHMT2D* to the *PYTHONPATH* environment variable. There are several ways to do it. See the user manual for details. One example is to add the path in your Python code.

.. code-block:: python

    import sys
    sys.path.append("/path/to/pyHMT2D")
..

Example Usage
-------------

There are at least two ways to use *pyHMT2D*:

- Command line interface (CLI): Only limited functions of *pyHMT2D* can be used in this way. You only need to
  type commands in a Windows terminal, e.g.,

  .. code-block:: bash

    $ hmt-calibrate calibration.json

  which runs a calibration job. Or

  .. code-block:: bash

    $ hmt-ras-to-srh Muncie2D.p01.hdf Terrain/TerrainMuncie_composite.tif srh_Muncie

  which converts a RAS 2D case to SRH-2D. See "examples/command_line_interface" for more details.

- Use in your own Python code (more flexibility)

  To use *pyHMT2D* in your Python code, simply add

  .. code-block:: python

   import pyHMT2D
  ..

  One example to use *pyHMT2D* to control the run of SRH-2D is as follows:

  .. code-block:: python

    #the follow should be modified based on your installation of SRH-2D
    version = "3.3"
    srh_pre_path = r"C:\Program Files\SMS 13.1 64-bit\Python36\Lib\site-packages\srh2d_exe\SRH_Pre_Console.exe"
    srh_path = r"C:\Program Files\SMS 13.1 64-bit\Python36\Lib\site-packages\srh2d_exe\SRH-2D_V330_Console.exe"
    extra_dll_path = r"C:\Program Files\SMS 13.1 64-bit\Python36\Lib\site-packages\srh2d_exe"

    #create a SRH-2D model instance
    my_srh_2d_model = pyHMT2D.SRH_2D.SRH_2D_Model(version, srh_pre_path,
                       srh_path, extra_dll_path, faceless=False)

    #initialize the SRH-2D model
    my_srh_2d_model.init_model()

    print("Hydraulic model name: ", my_srh_2d_model.getName())
    print("Hydraulic model version: ", my_srh_2d_model.getVersion())

    #open a SRH-2D project
    my_srh_2d_model.open_project("Muncie.srhhydro")

    #run SRH-2D Pre to preprocess the case
    my_srh_2d_model.run_pre_model()

    #run the SRH-2D model's current project
    my_srh_2d_model.run_model()

    #close the SRH-2D project
    my_srh_2d_model.close_project()

    #quit SRH-2D
    my_srh_2d_model.exit_model()
  ..

  Another example to use *pyHMT2D* to control the run of HEC-RAS is as follows:

  .. code-block:: python

    #create a HEC-RAS model instance
    my_hec_ras_model = pyHMT2D.RAS_2D.HEC_RAS_Model(version="5.0.7",
                                                    faceless=False)

    #initialize the HEC-RAS model
    my_hec_ras_model.init_model()

    #open a HEC-RAS project
    my_hec_ras_model.open_project("Muncie2D.prj",
                    "Terrain/TerrainMuncie_composite.tif")

    #run the HEC-RAS model's current project
    my_hec_ras_model.run_model()

    #close the HEC-RAS project
    my_hec_ras_model.close_project()

    #quit HEC-RAS
    my_hec_ras_model.exit_model()
  ..

  The last example is to use *pyHMT2D* to perform auto-calibration in two lines:

  .. code-block:: python

    my_calibrator = pyHMT2D.Calibration.Calibrator("calibration.json")

    my_calibrator.calibrate()
  ..

  More examples can be found in the "examples" directory.


Limitations
-----------

For SRH-2D:

-  This package is developed and tested with SRH-2D v3.3; other versions
   may work but has not been tested.
-  Currently, only flow data is processed; others such as sediment and water quality are ignored.
-  Currently *pyHMT2D* cannot manipulate other things such as hydraulic structures in the case configuration files.
   More functionalities will be added in the future.

For HEC-RAS 2D:

-  Only one 2D flow area is supported.
-  Only 2D flow area information is processed; others such as 1D
   channels and structures are ignored.
-  Currently, only flow data is processes; others such as sediment and water
   quality are ignored.
-  This package is developed and tested with HEC-RAS v5.0.7; other versions may work but have not been tested.

User manual and API documentation
---------------------------------

The *pyHMT2D* User Manual can be found in *docs*: `pyHMT2D_User_Manual.pdf <https://github.com/psu-efd/pyHMT2D/raw/main/docs/pyHMT2D_User_Manual.pdf>`__

The API documentation is hosted at `https://psu-efd.github.io/pyHMT2D_API_Web/ <https://psu-efd.github.io/pyHMT2D_API_Web/>`__

Acknowledgements and references
-------------------------------

*pyHMT2D* utilizes and/or benefits from several open source codes. The usage
of these codes strictly follows proper copyright laws and their licenses (see
the copies of their original licenses in the *licenses* directory). We
acknowledge their contributions.

In particular, the following packages were used and/or referenced:

-  `PyRAS - Python for River
   AnalysiS <https://github.com/solomonvimal/pyras>`__
-  `HaD-to-Py <https://github.com/latomkovic/HaD-to-Py>`__

Some of the examples and tests use dataset from public domain or authorized sources:

- Munice case data from HEC-RAS example data set (public domain)
- `Duck Pond <https://www.google.com/maps/@40.8041236,-77.8438126,522m/data=!3m1!1e3>`__ case data from Penn State University (with authorization for research and teaching purposes only)
- `Lidar data set from USGS <https://www.usgs.gov/core-science-systems/ngp/3dep>`_ (public domain)

The inclusion of these data sets in *pyHMT2D* is strictly for demonstration purpose only. Reuse or
repurpose of these dataset without explicit authorization from the original owner or copyright
holder is not permitted.

License
-------

MIT

Author
------

| Xiaofeng Liu, Ph.D., P.E.
| Associate Professor

| Department of Civil and Environmental Engineering
| Institute of Computational and Data Sciences
| Penn State University

223B Sackett Building, University Park, PA 16802

Web: http://water.engr.psu.edu/liu/

Contributors and contributor agreement
--------------------------------------
The list of contributors:
^^^^^^^^^^^^^^^^^^^^^^^^^
- (To be added)

Contributor agreement
^^^^^^^^^^^^^^^^^^^^^
First of all, thanks for your interest in contributing to *pyHMT2D*. Collectively, we can make *pyHMT2D* more
powerful, better, and easier to use.

Because of legal reasons and like many successful open source projects, contributors have to sign
a "Contributor License Agreement" to grant their rights to "Us". See details of the agreement on GitHub.
The signing of the agreement is automatic when a pull request is issued.

If you are just a user of *pyHMT2D*, the contributor agreement is irrelevant.
