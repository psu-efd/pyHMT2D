import numpy as np
import csv

import json

import pyHMT2D

from ..__common__ import pyHMT2D_SCALAR, pyHMT2D_VECTOR

class Calibrator(object):
    """Calibrator class to handle calibration process

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

        # configuration dictionary loaded from JSON
        self.configuration = {}

        self.load_config_json()

    def load_config_json(self):
        """ Load the JSON configuration file

        Returns
        -------

        """

        print("Load calibration configuration from file", self.config_json_file_name)

        with open(self.config_json_file_name) as f:
            self.configuration = json.load(f)

        #print(self.configuration)

        #print some information about the calibration configuration
        print("model =", self.configuration["model"])

        #model specific configuration
        print("Configuration for",self.configuration["model"], ":")
        print(json.dumps(self.configuration[self.configuration["model"]], indent=4, sort_keys=False))

        #objectives
        print("There are total", len(self.configuration["objectives"]), "calibration objectives:")
        print(json.dumps(self.configuration["objectives"], indent=4, sort_keys=False))