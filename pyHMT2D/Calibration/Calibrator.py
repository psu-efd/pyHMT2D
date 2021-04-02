import numpy as np

import json

import pyHMT2D

from ..__common__ import gVerbose

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

        # create a Backwater_1D_Model object
        if self.model_name == "Backwater-1D":
            self.hydraulic_model = pyHMT2D.Hydraulic_Models_Data.Backwater_1D_Model()

            # load Backwater_1D_Data configuration data
            self.hydraulic_data = pyHMT2D.Hydraulic_Models_Data.Backwater_1D_Data(self.config_json_file_name)

            # set the simulation case in the Backwater_1D_Model object
            self.hydraulic_model.set_simulation_case(self.hydraulic_data)


    def create_optimizer(self):
        """ Create the optimizer

        Returns
        -------

        """

        if self.optimizer_name == "scipy.optimize.local":
            self.optimizer = pyHMT2D.Calibration.Optimizer_ScipyOptimizeLocal(self.configuration["calibration"]["scipy.optimize.local"])

    def calibrate(self):
        """ Calibrate the model

        Returns
        -------

        """

        print("Calibration starts ...")

        if self.optimizer_name == "scipy.optimize.local":
            import scipy.optimize as OP


            # get the Manning's n calibration information
            materialID_list, initial_guess_list, ManningN_min_list, ManningN_max_list = \
                self.calibration_parameters.get_ManningN_Info_list()

            # build bounds
            ManningN_bounds = OP.Bounds(ManningN_min_list, ManningN_max_list)

            # scipy.optimize.minimize's arguments are different for different categories of methods. Need to
            # separate them.

            #for methods bounds are not allowed:
            if self.optimizer.method == "Nelder-Mead":
                result = OP.minimize(self.func_to_minimize, initial_guess_list, args=(materialID_list,),
                                 method=self.optimizer.method,
                                 tol=self.optimizer.tol,
                                 callback=self.callback,
                                 options=self.optimizer.options
                                 )
            #for methods bounds are allowd:
            elif self.optimizer.method == "L-BFGS-B":
                result = OP.minimize(self.func_to_minimize, initial_guess_list, args=(materialID_list,),
                                 method=self.optimizer.method,
                                 jac=self.optimizer.jac,
                                 hess=self.optimizer.hess,
                                 tol=self.optimizer.tol,
                                 bounds=ManningN_bounds,
                                 options=self.optimizer.options
                                 )

            print(result)

            #write out the optimizer's calibration intermediate steps' parameter values and calibration errors
            self.optimizer.write_optimization_results_to_csv()

            #write out the simulation results at measurement locations (for postprocessing)
            self.objectives.outputSimulationResultToCSV()

            print("Calibration ended.")

    def func_to_minimize(self, ManningNs, ManningN_MaterialIDs):
        """Function to minimize (the score, i.e., the cost function)

        Returns
        -------

        """

        total_score = np.inf

        if self.model_name == "Backwater-1D":
            #set the Manning's n with the new values
            for zoneI in range(len(ManningN_MaterialIDs)):
                self.hydraulic_model.get_simulation_case().modify_ManningsN(ManningN_MaterialIDs[zoneI],
                                                                          ManningNs[zoneI])

            # run the Backwater_1D_Model model
            self.hydraulic_model.run_model()

            # output the result to VTK
            vtkFileName = self.hydraulic_model.get_simulation_case().outputResultToVTK()

            #calculate the total score of the calibration run for all specified objectives
            #The score is calculated by comparing sampled result on VTK and measurement
            self.objectives.calculate_total_score(vtkFileName)

            if gVerbose: print("Total score = ", self.objectives.total_score)

            total_score = self.objectives.total_score

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

        "callback" is called by the optimizer every iteration, not every function evaluation call.

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
