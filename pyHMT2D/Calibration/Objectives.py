import numpy as np
import csv

import pyHMT2D

from ..__common__ import pyHMT2D_SCALAR, pyHMT2D_VECTOR

class Objective(object):
    """ Calibration objective class

    One objective corresponds to only one measurement (PointMeasurement, etc.). Its
    functionality is to sample on a 2D hydraulic modeling solution on the same measurement
    points or lines and calculate the objective score (norm of difference between
    measurement and simulation).

    Attributes
    ----------
        norm_order : str, optional
            order of norm to calculate the distance between measurement and simulation result (default = 2)

    """

    def __init__(self, measurement, norm_order=2):
        """Objective class constructor

        Parameters
        ----------
        measurement : Measurement object
            Object of the measurement (PointMeasurement, LineMeasurement, etc.)
        norm_order : str
            order of the norm to calculate the error
        """

        # Measurement object
        self.measurement = measurement

        # Order of the norm to calculate the error
        self.norm_order = norm_order

        # Objective score = norm of error
        self.score = 0.0

    def sample_on_result(self, vtkUnstructuredGridReader, varName):
        """Sample on result and calculate the objective score (norm of error)

        The sampling points are the same as in the measurement

        Parameters
        ----------
        vtkUnstructuredGridReader : vtkUnstructuredGridReader
            vtkUnstructuredGridReader object to pass along simulation result
        varName : str
            name of the variable to be sampled

        Returns
        -------

        """

        # Get the sampling points as vtkPoints
        sampling_points = self.measurement.get_measurement_points_as_vtkPoints()

        vtk_handler = pyHMT2D.Misc.vtkHandler()

        # sample on the sampling points
        points, varValues, elev_srh_2d = vtk_handler.probeUnstructuredGridVTKOverLine(
                                            sampling_points, vtkUnstructuredGridReader, varName)

        # calculate the difference
        if self.measurement.data_type == pyHMT2D_SCALAR:
            error = varValues - self.measurement.get_measurement_data()

            self.score = np.linalg.norm(error, self.norm_order)




class Objectives(object):
    """ Calibration objectives class

    An "Objectives" object is a list of "objective" objects.

    Attributes
    ----------

    """

    def __init__(self, name=""):
        """Objectives class constructor

        Parameters
        ----------
        name : str
            name of the Objectives object
        """

        # name of the Objectives
        self.name = name

        # list of all Objective objects
        self.objective_list = []
        