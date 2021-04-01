import numpy as np
import csv

import pyHMT2D

from ..__common__ import pyHMT2D_SCALAR, pyHMT2D_VECTOR

class Parameter(object):
    """ Calibration parameter base class


    Attributes
    ----------

    """

    def __init__(self):
        """Parameter class constructor

        Parameters
        ----------

        """

        self.type = "Parameter_Base"

class Parameter_ManningN(Parameter):
    """ ManningN calibration parameter class

    Attributes
    ----------
    materialID : int
        materialID for the ManningN
    initial_guess: float
        initial guess of the Manning's n
    min, max: float
        min and max of the calibration parameter bound

    """

    def __init__(self, parameterDict):
        """Parameter_ManningN class constructor

        Parameters
        ----------
        parameterDict : dict
            dictionary that contains ManningN calibration parameter information.

        """

        Parameter.__init__(self)

        # type
        self.type = "ManningN"

        # materialID of the Manning's n
        self.materialID = parameterDict["materialID"]

        # initial guess for the optimizer
        self.initial_guess = parameterDict["initial_guess"]

        # bound: min and max
        self.min = parameterDict["min"]
        self.max = parameterDict["max"]

        # some sanity check:
        if self.min > self.max or \
           self.initial_guess < self.min or \
           self.initial_guess > self.max:
            raise Exception("In ManningN calibration parameter specification, either the min is larger than"
                            " the max, or the initial guess is not in the [min, max] range. Check. Exiting ...")



class Parameters(object):
    """ Calibration parameter collection class

    A "Parameters" object is a list of "Parameter" objects.

    Attributes
    ----------

    """

    def __init__(self, parametersDict):
        """Objectives class constructor

        Parameters
        ----------
        parametersDict : dict
            a dictionary contains the information about calibration parameters

        """

        # dictionary containing information about calibration parameters
        self.parametersDict = parametersDict

        # list of all Parameter objects
        self.parameter_list = []

        # build parameter_list
        self.build_parameter_list()


    def build_parameter_list(self):
        """ Build parameter_list

        Returns
        -------

        """

        # loop through every calibration parameter in the dictionary
        for parameterDict in self.parametersDict:
            if parameterDict["type"] == "ManningN":
                #construct the Parameter object
                currParameter = pyHMT2D.Calibration.Parameter_ManningN(parameterDict)

                #append the Parameter_ManningN object to the list
                self.parameter_list.append(currParameter)
            else:
                raise Exception("The specified calibration parameter type", parameterDict["type"],
                                "is currently not supported. Support type is ManningN. Exiting ...")

    def get_ManningN_Info_list(self):
        """ Get the calibration ManningN's information

        This function returns the following:
        1. a list of calibration materialIDs
        2. a list of corresponding Manning's n initial guess values
        3. a list of corresponding Manning's n minimum value
        4. a list of corresponding Manning's n maximum value

        Returns
        -------
        materialID_list : list
            a list of calibration materialIDs
        initial_guess_list : list
            a list of Manning's n initial guesses
        ManningN_min_list : list
            a list of Manning's n minimum values
        ManningN_max_list : list
            a list of Manning's n maximum values


        """

        materialID_list = []
        initial_guess_list = []
        ManningN_min_list = []
        ManningN_max_list =[]


        #loop through all calibration parameter objects (whose type is "ManningN")
        for parameterI in self.parameter_list:
            if parameterI.type == "ManningN":
                materialID_list.append(parameterI.materialID)
                initial_guess_list.append(parameterI.initial_guess)
                ManningN_min_list.append(parameterI.min)
                ManningN_max_list.append(parameterI.max)

        #check the uniqueness of elements in materialID_list. If there are
        #repeated materialIDs, something is not right.
        if len(materialID_list) > len(set(materialID_list)):
            raise Exception("Manning's n calibration materailIDs are not unique. Please check. Exiting ...")

        return materialID_list, initial_guess_list, ManningN_min_list, ManningN_max_list
