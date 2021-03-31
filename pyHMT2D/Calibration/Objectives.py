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
        weight : float
            weight associated with this point measurement, in [0, 1]
        norm_order : str
            order of norm to calculate the distance between measurement and simulation result (default = 2)
        score : float
            score of the objective (norm of the error). The smaller the score, the closer the simulation
            result is to the measurement.

    """

    def __init__(self, solVarName, measurement, weight, norm_order=2):
        """Objective class constructor

        Parameters
        ----------
        solVarName : str
            Name of solution variable (used to sample the solution VTK file)
        measurement : Measurement object
            Object of the measurement (PointMeasurement, LineMeasurement, etc.)
        weight : float
            weight associated with this point measurement, in [0, 1]
        norm_order : str
            order of the norm to calculate the error
        """

        #solution variable (used to sample the solution VTK file)
        self.solVarName = solVarName

        # Measurement object
        self.measurement = measurement

        #check weight range in [0, 1]
        if weight < 0.0 or weight > 1.0:
            raise ValueError("Weight is not in the range of 0 and 1. Exiting...")

        self.weight = weight

        # Order of the norm to calculate the error
        self.norm_order = norm_order

        # Objective score = norm of error
        self.score = 0.0

    def calulate_score(self, vtkUnstructuredGridReader):
        """Sample on result and calculate the objective score (norm of error)

        The sampling points are the same as in the measurement

        Parameters
        ----------
        vtkUnstructuredGridReader : vtkUnstructuredGridReader
            vtkUnstructuredGridReader object to pass along simulation result

        Returns
        -------

        """

        # Get the sampling points as vtkPoints
        sampling_points = self.measurement.get_measurement_points_as_vtkPoints()

        vtk_handler = pyHMT2D.Misc.vtkHandler()

        # sample on the sampling points
        points, varValues, elev_srh_2d = vtk_handler.probeUnstructuredGridVTKOverLine(
                                            sampling_points, vtkUnstructuredGridReader,
                                            self.solVarName)

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

    def __init__(self, objectivesDict, name=""):
        """Objectives class constructor

        Parameters
        ----------
        objectivesDict : dict
            a dictionary contains the information about Objectives
        name : str, optional
            name of the Objectives object
        """

        # name of the Objectives
        self.name = name

        # dictionary containing information about Objectives
        self.objectivesDict = objectivesDict

        # total score (initialized as "inf"
        self.total_score = float('inf')

        # list of all Objective objects
        self.objective_list = []

        # build objective_list
        self.build_objective_list()


    def build_objective_list(self):
        """ Build objective_list

        Returns
        -------

        """

        # loop through every objective in the dictionary
        for objectiveDict in self.objectivesDict:
            #create a measurement associated with the current objective
            if objectiveDict["type"] == "PointMeasurement":

                #construct the PointMeasurement object
                currMeasurement = pyHMT2D.Calibration.PointMeasurement(objectiveDict["name"],
                                                                       objectiveDict["file"])

                #construct the Objective object
                currObjective = Objective(objectiveDict["solVarName"], currMeasurement, objectiveDict["weight"])

                #append the Objective object to the list
                self.objective_list.append(currObjective)


    def calculate_total_score(self, vtkUnstructuredGridFileName):
        """ Calculate the total score

        The total score is the weigthed summation of scores from all objectives

        Parameters
        ----------
        vtkUnstructuredGridFileName : str
            name of the vtkUnstructuredGrid file to pass along simulation result

        Returns
        -------

        """

        vtk_handler = pyHMT2D.Misc.vtkHandler()

        vtkUnstructuredGridReader = vtk_handler.readVTK_UnstructuredGrid(vtkUnstructuredGridFileName)

        self.total_score = 0.0

        total_weight = 0.0

        for objectiveI in self.objective_list:

            objectiveI.calulate_score(vtkUnstructuredGridReader)

            total_weight += objectiveI.weight

            self.total_score += objectiveI.score * objectiveI.weight

        self.total_score /= total_weight
