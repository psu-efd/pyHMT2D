import numpy as np
import csv

import json

import pyHMT2D

from ..__common__ import pyHMT2D_SCALAR, pyHMT2D_VECTOR

class Calibrator(object):
    """Calibrator class to handle calibration process

    A calibrator is constructed from its configuration file in JSON format. The configuration
    specifies hydraulic model (e.g., Backwater-1D, SRH-2D, HEC-RAS, etc.), case (model dependent), an Objectives
    object with a list of Objective objects, and optimizer.


    """

    def __init__(self, config_json_file_name):
        """ Calibrator class constructor

        Parameters
        ----------
        config_json_file_name : str
            name of the JSON configuration file for calibrator
        """

        # file name for the JSON configuration file
        self.config_json_file_name = config_json_file_name

        # model name, e.g., "SRH-2D", "HEC-RAS"
        self.model_name = ''

        # optimizer name, e.g., "scipy.optimize"
        self.optimizer_name = ''

        # configuration dictionary loaded from JSON
        self.configuration = {}

        # load the configuration from JSON
        self.load_config_json()

        # hydraulic model
        self.hydraulic_model = None

        # hydraulic data
        self.hydraulic_data = None

        #create hydraulic model and data
        self.create_model_data()

        # build the Objectives object
        self.objectives = pyHMT2D.Calibration.Objectives(self.configuration["objectives"])

    def load_config_json(self):
        """ Load the JSON configuration file

        Returns
        -------

        """

        print("Load calibration configuration from file", self.config_json_file_name)

        with open(self.config_json_file_name) as f:
            self.configuration = json.load(f)

        self.model_name = self.configuration["model"]

        self.optimizer_name = self.configuration["optimizer"]

        #print(self.configuration)

        #print some information about the calibration configuration
        print("model =", self.configuration["model"])

        #model specific configuration
        print("Configuration for",self.configuration["model"], ":")
        print(json.dumps(self.configuration[self.configuration["model"]], indent=4, sort_keys=False))

        #objectives
        print("There are total", len(self.configuration["objectives"]), "calibration objectives:")
        print(json.dumps(self.configuration["objectives"], indent=4, sort_keys=False))

    def create_model_data(self):
        """ Create the model and set up the case in preparation of calibration

        Returns
        -------

        """

        # create a Backwater_1D_Model object
        if self.model_name == "Backwater-1D":
            self.hydraulic_model = pyHMT2D.Hydraulic_Models_Data.Backwater_1D_Model()

            # load Backwater_1D_Data configuration data
            self.hydraulic_data = pyHMT2D.Hydraulic_Models_Data.Backwater_1D_Data(self.config_json_file_name)

            # set the simulation case in the Backwater_1D_Model object
            self.hydraulic_model.set_simulation_case(self.hydraulic_data)


    def calibrate(self):
        """ Calibrate the model

        Returns
        -------

        """

        if self.optimizer_name == "scipy.optimize":
            import scipy.optimize as OP


            #ordered Mateiral IDs corresonding to ManningN_initial
            ManningN_MaterialIDs = [1, 2]

            ManningN_initial = [0.02, 0.02]

            result = OP.minimize(self.func_to_minimize, ManningN_initial, args=(ManningN_MaterialIDs,), method='COBYLA',
                              options={'maxiter': 10000, 'rhobeg': 20})

    def func_to_minimize(self, ManningNs, ManningN_MaterialIDs):
        """Function to minimize (the score)

        Returns
        -------

        """

        if self.model_name == "Backwater-1D":
            #set the Manning's n with the new values
            for zoneI in range(len(ManningN_MaterialIDs)):
                self.hydraulic_model.get_simulation_case().modify_ManningsN(ManningN_MaterialIDs[zoneI],
                                                                          ManningNs[zoneI])

            # run the Backwater_1D_Model model
            self.hydraulic_model.run_model()

            # output the result to VTK
            vtkFileName = self.hydraulic_model.get_simulation_case().outputResultToVTK()

            #calculate the total score of the calibration run for all specified objectivies
            self.objectives.calculate_total_score(vtkFileName)

            print("Total score = ", self.objectives.total_score)

            return self.objectives.total_score





