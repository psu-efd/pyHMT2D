import numpy as np
from pathlib import Path
import os
import h5py

import json

import pyHMT2D

from ..__common__ import gVerbose

import logging

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

        # create hydraulic model and data
        self.create_model_data()

        # calibration parameters
        self.calibration_parameters = pyHMT2D.Calibration.Parameters(self.configuration["calibration"]["calibration_parameters"])

        # Objectives object
        self.objectives = pyHMT2D.Calibration.Objectives(self.configuration["calibration"]["objectives"])

        # Optimizer object
        self.optimizer = None

        # create the optimizer
        self.create_optimizer()

        # initialize a logger
        self.logger = None

        self.init_logger()


    def load_config_json(self):
        """ Load the JSON configuration file

        Returns
        -------

        """

        if gVerbose: print("Load calibration configuration from file", self.config_json_file_name)

        with open(self.config_json_file_name) as f:
            self.configuration = json.load(f)

        self.model_name = self.configuration["model"]

        self.optimizer_name = self.configuration["calibration"]["optimizer"]

        #print(self.configuration)

        #print some information about the calibration configuration
        if gVerbose: print("model =", self.configuration["model"])

        #model specific configuration
        if gVerbose: print("Configuration for",self.configuration["model"], ":")
        if gVerbose: print(json.dumps(self.configuration[self.configuration["model"]], indent=4, sort_keys=False))

        #calibration information
        #calibration parameters
        if gVerbose: print("There are total", len(self.configuration["calibration"]["calibration_parameters"]), "calibration parameters:")
        if gVerbose: print(json.dumps(self.configuration["calibration"]["calibration_parameters"], indent=4, sort_keys=False))

        #objectives
        if gVerbose: print("There are total", len(self.configuration["calibration"]["objectives"]), "calibration objectives:")
        if gVerbose: print(json.dumps(self.configuration["calibration"]["objectives"], indent=4, sort_keys=False))

        #optimizer
        if gVerbose: print("The selected optimizer is", self.configuration["calibration"]["optimizer"], "with the following setup:")
        if gVerbose: print(json.dumps(self.configuration["calibration"][self.configuration["calibration"]["optimizer"]], indent=4, sort_keys=False))


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

            # open the simulation case
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


    def create_optimizer(self):
        """ Create the optimizer

        Returns
        -------

        """

        if self.optimizer_name == "scipy.optimize.local":
            self.optimizer = pyHMT2D.Calibration.Optimizer_ScipyOptimizeLocal(self.configuration["calibration"]["scipy.optimize.local"])
        elif self.optimizer_name == "scipy.optimize.global":
            self.optimizer = pyHMT2D.Calibration.Optimizer_ScipyOptimizeGlobal(self.configuration["calibration"]["scipy.optimize.global"])
        elif self.optimizer_name == "enumerator":
            self.optimizer = pyHMT2D.Calibration.Optimizer_Enumerator(self.configuration["calibration"]["enumerator"])
        else:
            raise Exception("Specified optimizer is not supported.")

    def init_logger(self):
        """ Initialize a logger

        Returns
        -------

        """

        # Create a custom logger
        self.logger = logging.getLogger(__name__)

        # Create handlers
        c_handler = logging.StreamHandler()
        f_handler = logging.FileHandler('calibration.log', mode='w')

        # Create formatters and add it to handlers
        c_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
        f_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        c_handler.setFormatter(c_format)
        f_handler.setFormatter(f_format)

        # Add handlers to the logger
        self.logger.addHandler(c_handler)
        self.logger.addHandler(f_handler)

        self.logger.setLevel(logging.INFO)

    def calibrate(self):
        """ Calibrate the model

        Returns
        -------

        """

        print("Calibration starts ...")

        # get the Manning's n calibration information
        materialID_list, materialName_list, initial_guess_list, ManningN_min_list, ManningN_max_list = \
            self.calibration_parameters.get_ManningN_Info_list()

        if self.optimizer_name == "enumerator":

            self.optimizer.minimize(self.func_to_minimize, args=(materialID_list, materialName_list,), callback=self.callback)


        elif self.optimizer_name == "scipy.optimize.local":
            import scipy.optimize as OP

            # build bounds
            ManningN_bounds = OP.Bounds(np.array(ManningN_min_list), np.array(ManningN_max_list))

            # scipy.optimize.minimize's arguments are different for different categories of methods. Need to
            # separate them.

            #for methods bounds are not allowed:
            if self.optimizer.method == "Nelder-Mead":
                result = OP.minimize(self.func_to_minimize, initial_guess_list, args=(materialID_list,
                                                                                      materialName_list,),
                                 method=self.optimizer.method,
                                 callback=self.callback,
                                 options=self.optimizer.options
                                 )
            #for methods bounds are allowed:
            elif self.optimizer.method == "L-BFGS-B":
                result = OP.minimize(self.func_to_minimize, initial_guess_list, args=(materialID_list,
                                                                                      materialName_list,),
                                 method=self.optimizer.method,
                                 jac=self.optimizer.jac,
                                 hess=self.optimizer.hess,
                                 tol=self.optimizer.tol,
                                 bounds=ManningN_bounds,
                                 options=self.optimizer.options
                                 )
            elif self.optimizer.method == "Powell":
                result = OP.minimize(self.func_to_minimize, initial_guess_list, args=(materialID_list,
                                                                                      materialName_list,),
                                 method=self.optimizer.method,
                                 callback=self.callback,
                                 bounds=ManningN_bounds,
                                 options=self.optimizer.options
                                 )
            else:
                raise Exception("The optimization method %s is not supported." % self.optimizer.method)

            print(result)

        elif self.optimizer_name == "scipy.optimize.global":
            import scipy.optimize as OP

            # build ranges (tuples)
            ManningN_ranges = tuple(zip(ManningN_min_list, ManningN_max_list))

            #get the "full_output" option
            if self.configuration["calibration"]["scipy.optimize.global"]["full_output"] == "True":
                full_output = True
            elif self.configuration["calibration"]["scipy.optimize.global"]["full_output"] == "False":
                full_output = False
            else:
                raise Exception("Optinmization parameter full_output can only be True or False. Please check.")

            if self.optimizer.method == "brute":
                result = OP.brute(self.func_to_minimize, ManningN_ranges, args=(materialID_list,
                                                                                materialName_list,),
                                  Ns=self.configuration["calibration"]["scipy.optimize.global"]["Ns"],
                                  #finish=self.configuration["calibration"]["scipy.optimize.global"]["finish"],
                                  finish=None,   #hard-wired here to not allow "brute" to go beyond the range
                                  full_output=full_output
                                 )

                # send information to logger for the record or restart
                if full_output: #return both parameter values and function value
                    msg = "\nThe optimized parameters:: " + ", ".join(map(str, result[0])) \
                          + "\nThe function value at optimized parameters:  " + str(result[1])
                else:          #return only paramter values
                    msg = "\nThe optimized parameters:: " + ", ".join(map(str, result))

                self.logger.info(msg)

            else:
                raise Exception("The optimization method %s is not supported." % self.optimizer.method)

        else:
            raise  Exception("The optimizer name %s is not supported." % self.optimizer_name)



        #write out the optimizer's calibration intermediate steps' parameter values and calibration errors
        self.optimizer.write_optimization_results_to_csv()

        #write out the simulation results at measurement locations (for postprocessing)
        self.objectives.outputSimulationResultToCSV()

        print("Calibration ended.")

    def func_to_minimize(self, ManningNs, ManningN_MaterialIDs, ManningN_MaterialNames):
        """Function to minimize (the score, i.e., the cost function)

        Returns
        -------

        """

        total_score = np.inf

        if self.model_name == "Backwater-1D":
            #set the Manning's n with the new values
            self.hydraulic_model.get_simulation_case().modify_ManningsN(ManningN_MaterialIDs,
                                                                        ManningNs,
                                                                        ManningN_MaterialNames)

            # run the Backwater_1D_Model model
            self.hydraulic_model.run_model()

            # output the result to VTK
            vtkFileName = self.hydraulic_model.get_simulation_case().outputResultToVTK()

            #calculate the total score of the calibration run for all specified objectives
            #The score is calculated by comparing sampled result on VTK and measurement
            self.objectives.calculate_total_score(vtkFileName)

            if gVerbose: print("Total score = ", self.objectives.total_score)

            total_score = self.objectives.total_score

        elif self.model_name == "SRH-2D":
            # set the Manning's n with the new values
            self.hydraulic_data.srhhydro_obj.modify_ManningsN(ManningN_MaterialIDs,
                                                              ManningNs,
                                                              ManningN_MaterialNames)

            srhhydro_filename = self.hydraulic_data.srhhydro_obj.srhhydro_filename

            self.hydraulic_data.srhhydro_obj.write_to_file(srhhydro_filename)

            # run SRH-2D Pre to preprocess the case
            self.hydraulic_model.run_pre_model()

            # run the SRH-2D model's current project
            self.hydraulic_model.run_model()

            # output the result to VTK
            # read SRH-2D result in XMDF format (*.h5)
            # Whether the XMDF result is nodal or cell center. In SRH-2D's ".srhhydro" file,
            # the output option for "OutputFormat" can be manually changed before simulation.
            # Options are "XMDF" (results at at nodes), "XMDFC" (results are at cell centers), etc.
            # For example, "OutputFormat XMDFC EN". The following lines show that the SRH-2D simulation
            # was run with "XMDFC" as output format (see the "XMDFC" part of the result file name) and thus
            # we set "bNodal = False".
            bNodal = False

            self.hydraulic_data.readSRHXMDFFile(self.hydraulic_data.get_case_name() + "_XMDFC.h5", bNodal)

            # export the SRH-2D result to VTK: lastTimeStep=True means we only want to deal with the last time step.
            # See the code documentation of outputXMDFDataToVTK(...) for more options. It returns a list of vtk file names.
            vtkFileNameList = self.hydraulic_data.outputXMDFDataToVTK(bNodal, lastTimeStep=True, dir='')

            # calculate the total score of the calibration run for all specified objectives
            # The score is calculated by comparing sampled result on VTK and measurement
            self.objectives.calculate_total_score(vtkFileNameList[-1])  #only take the last vtk result file

            if gVerbose: print("Total score = ", self.objectives.total_score)

            total_score = self.objectives.total_score

        elif self.model_name == "HEC-RAS":
            # open the simulation case
            self.hydraulic_model.open_project(self.configuration["HEC-RAS"]["case"],
                                              self.configuration["HEC-RAS"]["terrainFileName"])

            self.hydraulic_data = self.hydraulic_model.get_simulation_case()

            # set the Manning's n with the new values
            self.hydraulic_data.modify_ManningsN(ManningN_MaterialIDs,
                                                            ManningNs,
                                               ManningN_MaterialNames)

            #update the time stamp of the Manning's n GeoTiff file (to force HEC-RAS to re-compute 2D flow area's
            #properties table. (No need? The above Manning's modification already updatet the time stamp.)
            if os.path.dirname(self.hydraulic_data.hdf_filename) == '':
                fileBase = b''
            else:
                fileBase = str.encode(os.path.dirname(self.hydraulic_data.hdf_filename) + '/')

            full_landcover_filename = (fileBase + self.hydraulic_data.landcover_filename).decode("ASCII")

            Path(full_landcover_filename).touch()

            #save the current project before run it
            self.hydraulic_model.save_project()

            # run the HEC-RAS model's current project
            self.hydraulic_model.run_model()

            # read the HEC-RAS simulation result
            self.hydraulic_data.load2DAreaSolutions()

            # save the HEC-RAS simulation result to VTK. It returns a list of VTK files
            vtkFileNameList = self.hydraulic_data.saveHEC_RAS2D_results_to_VTK(lastTimeStep=True)

            # calculate the total score of the calibration run for all specified objectives
            # The score is calculated by comparing sampled result on VTK and measurement
            self.objectives.calculate_total_score(vtkFileNameList[-1])  #only take the last vtk result file

            if gVerbose: print("Total score = ", self.objectives.total_score)

            total_score = self.objectives.total_score

            #close project
            self.hydraulic_model.close_project()


        else:
            raise Exception("The specified model: %s, is not supported", self.model_name)


        #send information to logger for the record or restart
        msg =   "ManningN_MaterialIDs: " + ", ".join(map(str, ManningN_MaterialIDs))\
              + " ManningNs: " + ", ".join(map(str, ManningNs)) \
              + " total_score: " + str(total_score)
        self.logger.info(msg)


        # updating the optimizer's lists for record. Pass to optimizer without arguments or parentheses.
        # ref: https://stackoverflow.com/questions/16739065/how-to-display-progress-of-scipy-optimize-function
        if not self.optimizer.num_calls:  # first call is stored in all lists
            self.optimizer.decreasing_list_calls_inp.append(ManningNs)
            self.optimizer.decreasing_list_calls_res.append(total_score)
            self.optimizer.list_callback_inp.append(ManningNs)
            self.optimizer.list_callback_res.append(total_score)
        elif total_score < self.optimizer.decreasing_list_calls_res[-1]:
            self.optimizer.decreasing_list_calls_inp.append(ManningNs)
            self.optimizer.decreasing_list_calls_res.append(total_score)

        self.optimizer.list_calls_inp.append(ManningNs)
        self.optimizer.list_calls_res.append(total_score)
        self.optimizer.num_calls += 1

        return total_score

    #ref: https://stackoverflow.com/questions/16739065/how-to-display-progress-of-scipy-optimize-function
    def callback(self, xk, *_):
        """Callback function that can be used by optimizers of scipy.optimize.

        "callback" is called by the optimizer every iteration, not every function evaluation call. For "enumerator"
        optimizer, it is called for each parameter combination.

        The third argument "*_" makes sure that it still works when the
        optimizer calls the callback function with more than one argument. Pass
        to optimizer without arguments or parentheses.

        """

        s1 = ""
        xk = np.atleast_1d(xk)

        # search backwards in input list for input corresponding to xk
        for i, x in reversed(list(enumerate(self.optimizer.list_calls_inp))):
            x = np.atleast_1d(x)
            if np.allclose(x, xk):
                break

        for comp in xk:
            s1 += f"{comp:10.5e}\t"
        s1 += f"{self.optimizer.list_calls_res[i]:10.5e}"

        self.optimizer.list_callback_inp.append(xk)
        self.optimizer.list_callback_res.append(self.optimizer.list_calls_res[i])

        # if first call, print out the header.
        if not self.optimizer.callback_count:
            s0 = ""
            for j, _ in enumerate(xk):
                tmp = f"Parameter-{j + 1}"
                s0 += f"{tmp:10s}\t"
            s0 += "Calibration-Error"
            print(s0)

        # print the current iteration parameter values and calibration error
        print(s1)

        self.optimizer.callback_count += 1

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
