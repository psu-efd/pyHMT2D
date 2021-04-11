.. image:: https://github.com/psu-efd/pyHMT2D/raw/main/docs/source/_static/images/pyhmt2d_logo_vel_color_vector.png
   :target: https://github.com/psu-efd/pyHMT2D
   :alt: pyHMT2D

*pyHMT2D* - Python Hydraulic Modeling Tools - 2D
================================================

*pyHMT2D* is a Python package developed to control and (semi)automate 2D
hydraulic modeling, and pre-/postprocessing simulation results.
Currently, the following hydraulic models are supported:

-  SRH-2D
-  HEC-RAS

In the future, suport for more 2D models will be added.

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

-  read SRH-2D results
-  convert SRH-2D results to VTK format, one of the most popular format for
   scientific data
-  sample and probe simulation results (with the functionalities of VTK
   library)
-  control SRH-2D simulations

For HEC-RAS 2D modeling:

-  read RAS2D results (HDF files)
-  convert RAS2D results to VTK

   -  point and cell center data (depht, water surface elevation,
      velocity, etc.)
   -  interpolate between point and cell data
   -  face data (e.g., subterrain data)

-  convert RAS2D mesh, boundary conditions, and Manningâ€™s n data into
   SRH-2D format such that SRH-2D and HEC-RAS 2D can run a case with
   exactly the same mesh for comparison purpose. Additionally, HEC-RAS
   2D can be used as a mesh generator for SRH-2D.
-  sample and probe simulation results (with the functionalities of VTK
   library)
-  control HEC-RAS 2D simulations

With the control and automation capability above, it is much easier to
do the following:

-  automatic calibration of models with available optimization and
   calibration Python packages. Currently, "scipy"'s "optimize" module is supported, which
   includes many local and global optimization methods.
-  Monte-Carlo simulations with scripting and Pythonâ€™s statistic
   libraries
-  ...

Other features:

-  calculate the difference between simulation results (regardless they are on the same mesh or not)
-  create and manipulate terrain data

Requirements
------------

This package uses the following libraries:

-  `h5py <https://www.h5py.org/>`__
-  `vtk <https://github.com/Kitware/VTK>`__
-  `pywin32 <https://pypi.org/project/pywin32/>`__
-  `gdal <https://pypi.org/project/GDAL/>`__
-  `affine <https://pypi.org/project/affine/>`__

You can install these requirements all at once by first getting the the "`requirements.txt <https://github.com/psu-efd/pyHMT2D/blob/main/requirements.txt>`__"
file from *pyHMT2D*'s GitHub Repository. Then, in a terminal (preferably with a dedicated Python virtual environment activated), do

.. code-block:: bash

   $ pip install -r requirements.txt

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

     If you use this approach to install *pyHMT2D*, in order to make the Python Interpreter aware of *pyHMT2D*, you need to add the path to *pyHMT2D* to the *PYTHONPATH*. There are several ways to do it. See the user manual for details. One example is to add the path in your Python code.

.. code-block:: python

    import sys
    sys.path.append("/path/to/pyHMT2D")
..

Example Usage
-------------

To use *pyHMT2D* in your Python code, simply add

.. code-block:: python

    import pyHMT2D
..

One example to use *pyHMT2D* to control the run of HEC-RAS is as follows:

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

More examples can be found in the "examples" directory.


Limitations
-----------

For SRH-2D:

-  This package is developed and tested with SRH-2D v3.3; other versions
   may work but has not been tested.
-  Currently *pyHMT2D* cannot manipulate other things such as hydraulic structures in the case configuration files.
   More functionalities will be added in the future.

For HEC-RAS 2D:

-  Only one 2D flow area is supported.
-  Only 2D flow area information is processed; others such as 1D
   channels and structures are ignored.
-  Only flow data is processes; others such as sediment and water
   quality are ignored.
-  This package is developed and tested with HEC-RAS v5.0.7 and v6.0
   beta; other versions may work but has not been tested.

Acknowledgements and references
-------------------------------

*pyHMT2D* utilizes and/or benefits from several open source codes. The usage
of these codes strictly follows proper copyright laws and their licenses (see
the copies of their original licenses in the â€œlicensesâ€ directory). We
acknowledge their contributions.

In particular, the following packages were used:

-  `PyRAS - Python for River
   AnalysiS <https://github.com/solomonvimal/pyras>`__
-  `HaD-to-Py <https://github.com/latomkovic/HaD-to-Py>`__

Some of the examples and tests use dataset from public domain or authorized sources:

- Munice case data from HEC-RAS example data set (public domain)
- Duck Pond case data from Penn State University (with authorization for research and teaching purposes only)
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
