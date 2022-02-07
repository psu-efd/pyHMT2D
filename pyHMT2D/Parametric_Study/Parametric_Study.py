"""
Parametric Study.

Sampling in the parameter space.

"""
import os
import json
import copy
import shutil

from pyHMT2D.Calibration.Parameters import *

from pyHMT2D.Parametric_Study.Sampler import *

from ..__common__ import gVerbose

import logging

class Parametric_Study(object):
    """ Parametric study base class


    Attributes
    ----------

    """

    def __init__(self, config_json_file_name):
        """Parametric study class constructor

        Parameters
        ----------
        config_json_file_name : str
            file name for the json configuration file

        """

        # file name for the JSON configuration file
        self.config_json_file_name = config_json_file_name

        # model name, e.g., "SRH-2D", "HEC-RAS"
        self.model_name = ''

        # sampler name, e.g., "skopt.sampler"
        self.sampler_name = ''

        # configuration dictionary loaded from JSON
        self.configuration = {}

        # load the configuration from JSON
        self.load_config_json()

        # hydraulic model
        self.hydraulic_model = None

        # hydraulic data
        self.hydraulic_data = None

        # create hydraulic model and data
        self.create_model_data()

        # parameters for the parametric study
        self.parameters = Parameters(self.configuration["parametric_study"]["parameters"])

        # number of samples
        self.n_samples = self.configuration["parametric_study"]["n_samples"]

        # sampler object
        self.sampler = None

        # create the sampler
        self.create_sampler()

        # initialize a logger
        self.logger = None

        self.init_logger()


    def load_config_json(self):
        """ Load the JSON configuration file

        Returns
        -------

        """

        if gVerbose: print("Load parametric study configuration from file", self.config_json_file_name)

        with open(self.config_json_file_name) as f:
            self.configuration = json.load(f)

        self.model_name = self.configuration["model"]

        self.sampler_name = self.configuration["parametric_study"]["sampler"]

        #print(self.configuration)

        #print some information about the calibration configuration
        if gVerbose: print("model =", self.configuration["model"])

        #model specific configuration
        if gVerbose: print("Configuration for",self.configuration["model"], ":")
        if gVerbose: print(json.dumps(self.configuration[self.configuration["model"]], indent=4, sort_keys=False))

        #parametric study information
        # parameters
        if gVerbose: print("There are total", len(self.configuration["parametric_study"]["parameters"]), " parameters:")
        if gVerbose: print(json.dumps(self.configuration["parametric_study"]["parameters"], indent=4, sort_keys=False))

        #sampler
        if gVerbose: print("The selected sampler is", self.configuration["parametric_study"]["sampler"], "with the following setup:")
        if gVerbose: print(json.dumps(self.configuration["parametric_study"][self.configuration["parametric_study"]["sampler"]], indent=4, sort_keys=False))


    def create_model_data(self):
        """ Create the model and set up the case in preparation of calibration

        Returns
        -------

        """

        if self.model_name == "Backwater-1D": # create a Backwater_1D_Model object and open the simulation case
            self.hydraulic_model = pyHMT2D.Hydraulic_Models_Data.Backwater_1D_Model()

            # load Backwater_1D_Data configuration data
            self.hydraulic_data = pyHMT2D.Hydraulic_Models_Data.Backwater_1D_Data(self.config_json_file_name)

            # set the simulation case in the Backwater_1D_Model object
            self.hydraulic_model.set_simulation_case(self.hydraulic_data)

        elif self.model_name == "SRH-2D": # create a SRH_2D_Model object and open the simulation case
            version = self.configuration["SRH-2D"]["version"]
            srh_pre_path = self.configuration["SRH-2D"]["srh_pre_path"]
            srh_path = self.configuration["SRH-2D"]["srh_path"]
            extra_dll_path = self.configuration["SRH-2D"]["extra_dll_path"]

            # create a SRH-2D model instance
            self.hydraulic_model = pyHMT2D.SRH_2D.SRH_2D_Model(version, srh_pre_path,
                                                          srh_path, extra_dll_path, faceless=False)

            # initialize the SRH-2D model
            self.hydraulic_model.init_model()

            if gVerbose: print("Hydraulic model name: ", self.hydraulic_model.getName())
            if gVerbose: print("Hydraulic model version: ", self.hydraulic_model.getVersion())

            # open the template case
            self.hydraulic_model.open_project(self.configuration["SRH-2D"]["case"])

            self.hydraulic_data = self.hydraulic_model.get_simulation_case()

        elif self.model_name == "HEC-RAS": # create a HEC_RAS_Model object and open the simulation case
            version = self.configuration["HEC-RAS"]["version"]

            # whether to run HEC-RAS faceless
            if (self.configuration["HEC-RAS"]["faceless"] == "True"):
                faceless = True
            elif (self.configuration["HEC-RAS"]["faceless"] == "False"):
                faceless = False
            else:
                raise Exception("faceless should be either True or False. Please check.")

            # create a HEC-RAS model instance
            self.hydraulic_model = pyHMT2D.RAS_2D.HEC_RAS_Model(version, faceless)

            # initialize the HEC-RAS model
            self.hydraulic_model.init_model()

            print("Hydraulic model name: ", self.hydraulic_model.getName())
            print("Hydraulic model version: ", self.hydraulic_model.getVersion())

            # open the simulation case
            #self.hydraulic_model.open_project(self.configuration["HEC-RAS"]["case"],
            #                                  self.configuration["HEC-RAS"]["terrainFileName"])

            #self.hydraulic_data = self.hydraulic_model.get_simulation_case()


        else:
            raise Exception("The specified model: %s, is not supported", self.model_name)


    def create_sampler(self):
        """ Create the sampler

        Returns
        -------

        """

        if self.sampler_name == "skopt.sampler":
            self.sampler = Sampler_SkoptSampler(self.configuration["parametric_study"]["skopt.sampler"])
        elif self.samplerr_name == "enumerator":
            self.sampler = Sampler_Enumerator(self.configuration["parametric_study"]["enumerator"])
        else:
            raise Exception("Specified sampler is not supported.")

    def init_logger(self):
        """ Initialize a logger

        Returns
        -------

        """

        # Create a custom logger
        self.logger = logging.getLogger(__name__)

        # Create handlers
        c_handler = logging.StreamHandler()
        f_handler = logging.FileHandler('parametric_study.log', mode='w')

        # Create formatters and add it to handlers
        c_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
        f_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        c_handler.setFormatter(c_format)
        f_handler.setFormatter(f_format)

        # Add handlers to the logger
        self.logger.addHandler(c_handler)
        self.logger.addHandler(f_handler)

        self.logger.setLevel(logging.INFO)

    def create_all_cases(self):
        """ Create all cases

        Returns
        -------

        """

        print("Create all cases ...")

        # all cases will be in the "cases" directory (make sure delete this directory before run the script)
        if os.path.isdir("cases"):
            raise Exception("The directory cases already exists. Make sure to remove it before run this program.")
        else:
            os.mkdir("cases")

        #build the parameter space
        param_space_list = self.parameters.get_parameter_space()
        param_space = Space(param_space_list)

        samples = self.sampler.generate(param_space, self.n_samples)
        print("Samples: ", samples)

        if self.model_name == "SRH-2D":

            #loop over all samples
            for i, sample in zip(range(len(samples)), samples):
                print("current sample: ", i, sample)

                sample_np = np.array(sample)

                my_srh_2d_srhhydro = copy.deepcopy(self.hydraulic_data.srhhydro_obj)

                srh_caseName = my_srh_2d_srhhydro.get_Case_Name() + "_" + str(i)

                # new grid file name, HydroMat file name, and srhhydro file name
                newGridFileName = srh_caseName + ".srhgeom"
                newHydroMatFileName = srh_caseName + ".srhmat"
                newSRHHydroFileName = srh_caseName + ".srhhydro"

                # old grid file name, HydroMat file name
                oldGridFileName = my_srh_2d_srhhydro.get_Grid_FileName()
                oldHydroMatFileName = my_srh_2d_srhhydro.get_HydroMat_FileName()

                #copy to old grid file and HydroMat files to new ones (no need to change anything there)
                shutil.copy(oldGridFileName, newGridFileName)
                shutil.copy(oldHydroMatFileName, newHydroMatFileName)

                my_srh_2d_srhhydro.modify_Case_Name(srh_caseName)
                my_srh_2d_srhhydro.modify_Grid_FileName(newGridFileName)
                my_srh_2d_srhhydro.modify_HydroMat_FileName(newHydroMatFileName)

                # go over each parameter in the current sample
                for parameterI in range(len(sample)):
                    if self.parameters.parameter_list[parameterI].type == "ManningN":  # need to take care of active or not later; now assume all active
                        materialIDs = [self.parameters.parameter_list[parameterI].materialID]
                        newManningsNValues = [sample_np[parameterI]]
                        ManningN_MaterialNames = [self.parameters.parameter_list[parameterI].name]

                        my_srh_2d_srhhydro.modify_ManningsN(materialIDs, newManningsNValues, ManningN_MaterialNames)

                    elif self.parameters.parameter_list[parameterI].type == "InletQ":
                        bcIDs = [self.parameters.parameter_list[parameterI].bcID]
                        newInletQValues = [sample_np[parameterI]]

                        my_srh_2d_srhhydro.modify_InletQ(bcIDs, newInletQValues)

                # save the srhhydro file
                my_srh_2d_srhhydro.write_to_file(newSRHHydroFileName)

                # make a case directory inside "cases"
                os.mkdir("cases/case_" + str(i))

                # move the three case files to the current case directory
                shutil.move(newSRHHydroFileName, "cases/case_" + str(i))
                shutil.move(newHydroMatFileName, "cases/case_" + str(i))
                shutil.move(newGridFileName, "cases/case_" + str(i))

        else:
            raise NotImplementedError


        return samples, self.parameters.get_parameter_name_list()


    def close_and_cleanup(self):
        """Close the calibration process and clean up

        Returns
        -------

        """

        if self.model_name == "Backwater-1D":
            pass
        elif self.model_name == "SRH-2D":
            pass
        elif self.model_name == "HEC-RAS":
            # quit HEC-RAS
            self.hydraulic_model.exit_model()
