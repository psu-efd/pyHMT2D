
from ..Misc import json_dict_type_correction

import numpy as np

class Optimizier(object):
    """ Optimizer base class


    Attributes
    ----------

    """

    def __init__(self, optimizerDict):
        """Optimizer class constructor

        Parameters
        ----------
        optimizerDict : dict
            dictionary that contains information on the optimizer.

        """

        self.type = "Optimizer_Base"

        self.optimizerDict = optimizerDict

    def load_from_optimizer_dict(self):
        """ Load configuration of the optimizer from dictionary

        Returns
        -------

        """

        pass

class Optimizer_ScipyOptimizeLocal(Optimizier):
    """ Optimizer using scipy.optimize.local class

    Using scipy.optimize's local optimization, i.e., use the "minimize(...)" function

    Attributes
    ----------

    """

    def __init__(self, optimizerDict):
        """Optimizer_ScipyOptimizeLocal class constructor

        Parameters
        ----------

        """

        Optimizier.__init__(self, optimizerDict)

        # type
        self.type = "ScipyOptimizeLocal"

        # method: default "L-BFGS-B". The specified method should support
        # the specification of "bounds".
        self.method = "L-BFGS-B"

        # Method for computing the gradient vector. Only for CG, BFGS, Newton-CG, L-BFGS-B, TNC,
        # SLSQP, dogleg, trust-ncg, trust-krylov, trust-exact and trust-constr.
        self.jac = None

        # Method for computing the Hessian matrix. Only for Newton-CG, dogleg, trust-ncg,
        # trust-krylov, trust-exact and trust-constr.
        self.hess = None

        # Maximum function evalutions (only for a subsset of methods, e.g., L-BFGS-B)
        self.maxfun = 1000

        # Tolerance for termination. For detailed control, use solver-specific options.
        self.tol = 0.001

        #options:

        # Maximum function iteraions (only for a subset of methods, e.g., L-BFGS-B)
        self.maxiter = 2000

        # Set to True to print convergence messages.
        self.disp = True

        # dictionary for the options of the particular method
        self.options = {}

        # load optimizer configuration from dictionary
        self.load_from_optimizer_dict()

        # for callback: record the data during the optimization such as number of function calls,
        # ref: https://stackoverflow.com/questions/16739065/how-to-display-progress-of-scipy-optimize-function
        self.num_calls = 0  # how many times the cost function has been called
        self.callback_count = 0  # number of times callback has been called (= iteration count)
        self.list_calls_inp = []  # input of all calls
        self.list_calls_res = []  # result of all calls
        self.decreasing_list_calls_inp = []  # input of calls that resulted in decrease
        self.decreasing_list_calls_res = []  # result of calls that resulted in decrease
        self.list_callback_inp = []  # only appends inputs on callback, as such they correspond to the iterations
        self.list_callback_res = []  # only appends results on callback, as such they correspond to the iterations

    def load_from_optimizer_dict(self):
        """ Load configuration of the optimizer from the dictionary

        Returns
        -------

        """

        if "method" in self.optimizerDict:
            self.method = self.optimizerDict["method"]

        if "jac" in self.optimizerDict:
            self.jac = self.optimizerDict["jac"]

        if "hess" in self.optimizerDict:
            self.hess = self.optimizerDict["hess"]

        if "maxfun" in self.optimizerDict:
            self.maxfun = self.optimizerDict["maxfun"]

        if "maxiter" in self.optimizerDict:
            self.maxiter = self.optimizerDict["maxiter"]

        if "tol" in self.optimizerDict:
            self.tol = self.optimizerDict["tol"]

        if "disp" in self.optimizerDict:
            if self.optimizerDict["disp"] == "True":
                self.disp = True
            elif self.optimizerDict["disp"] == "False":
                self.disp = False
            else:
                raise Exception("Optimizer's disp option should be either True of False. Please check. Exiting ...")

        if "options" in self.optimizerDict:
            self.options = self.optimizerDict["options"]

            #correct
            json_dict_type_correction(self.options)


    def write_optimization_results_to_csv(self, parameterNames=[]):
        """ Write optimization results to files in csv format

        Returns
        -------

        """

        #write out the calibration parameters' values and calibration errors for all iterations
        fid = open(self.type+"_"+self.method+"_calibration_results"+".csv", "w")

        #loop through all parameters
        for parameterI in range(len(self.list_callback_inp[0])):
            if len(parameterNames) == 0:
                fid.write("Parameter-%d," % parameterI)
            else:
                fid.write("%s," % parameterNames[parameterI])

        fid.write("Calibration-Error\n")

        #loop through all iterations
        for iterationI in range(len(self.list_callback_inp)):
            for parameterI in range(len(self.list_callback_inp[0])):
                fid.write("%10.5e," % self.list_callback_inp[iterationI][parameterI])

            fid.write("%10.5e\n" % self.list_callback_res[iterationI])

        fid.close()


class Optimizer_ScipyOptimizeGlobal(Optimizier):
    """ Optimizer using scipy.optimize.global class

    Using scipy.optimize's global optimization, e.g., use the "brute(...)" function

    Attributes
    ----------

    """

    def __init__(self, optimizerDict):
        """Optimizer_ScipyOptimizeGlobal class constructor

        Parameters
        ----------

        """

        Optimizier.__init__(self, optimizerDict)

        # type
        self.type = "ScipyOptimizeGlobal"

        # method: default "brute".
        self.method = "brute"

        #options:

        # dictionary for the options of the particular method
        self.options = {}

        # load optimizer configuration from dictionary
        self.load_from_optimizer_dict()

        # for callback: record the data during the optimization such as number of function calls,
        # ref: https://stackoverflow.com/questions/16739065/how-to-display-progress-of-scipy-optimize-function
        self.num_calls = 0  # how many times the cost function has been called
        self.callback_count = 0  # number of times callback has been called (= iteration count)
        self.list_calls_inp = []  # input of all calls
        self.list_calls_res = []  # result of all calls
        self.decreasing_list_calls_inp = []  # input of calls that resulted in decrease
        self.decreasing_list_calls_res = []  # result of calls that resulted in decrease
        self.list_callback_inp = []  # only appends inputs on callback, as such they correspond to the iterations
        self.list_callback_res = []  # only appends results on callback, as such they correspond to the iterations

    def load_from_optimizer_dict(self):
        """ Load configuration of the optimizer from the dictionary

        Returns
        -------

        """

        if "method" in self.optimizerDict:
            self.method = self.optimizerDict["method"]

        if "options" in self.optimizerDict:
            self.options = self.optimizerDict["options"]

            #correct
            json_dict_type_correction(self.options)


    def write_optimization_results_to_csv(self, parameterNames=[]):
        """ Write optimization results to files in csv format

        Returns
        -------

        """

        #write out the calibration parameters' values and calibration errors for all iterations
        fid = open(self.type+"_"+self.method+"_calibration_results"+".csv", "w")

        #loop through all parameters
        for parameterI in range(len(self.list_callback_inp[0])):
            if len(parameterNames) == 0:
                fid.write("Parameter-%d," % parameterI)
            else:
                fid.write("%s," % parameterNames[parameterI])

        fid.write("Calibration-Error\n")

        #loop through all iterations
        for iterationI in range(len(self.list_callback_inp)):
            for parameterI in range(len(self.list_callback_inp[0])):
                fid.write("%10.5e," % self.list_callback_inp[iterationI][parameterI])

            fid.write("%10.5e\n" % self.list_callback_res[iterationI])

        fid.close()

class Optimizer_Enumerator(Optimizier):
    """ Optimizer using user-provided enumeration of parameter combinations class

    User needs to provide a list of parameter combinations. Then, the algorithm sweeps
    through the enumeration and picks the best one.

    Attributes
    ----------

    """

    def __init__(self, optimizerDict):
        """Optimizer_Enumerator class constructor

        Parameters
        ----------

        """

        Optimizier.__init__(self, optimizerDict)

        # type
        self.type = "enumerator"

        # method: default (not relevant).
        self.method = "None"

        #options:

        # dictionary for the options of the particular method
        self.options = {}

        #parameter_combinations (a list of lists)
        #For each list, it contains: "materialID": Manning's n dictionary
        self.parameter_combinations = []

        # load optimizer configuration from dictionary
        self.load_from_optimizer_dict()

        #best combination index
        self.best_combination_index = -1

        #best calibration score
        self.best_calibration_score = np.inf

        # for callback: record the data during the optimization such as number of function calls,
        # ref: https://stackoverflow.com/questions/16739065/how-to-display-progress-of-scipy-optimize-function
        self.num_calls = 0  # how many times the cost function has been called
        self.callback_count = 0  # number of times callback has been called (= iteration count)
        self.list_calls_inp = []  # input of all calls
        self.list_calls_res = []  # result of all calls
        self.decreasing_list_calls_inp = []  # input of calls that resulted in decrease
        self.decreasing_list_calls_res = []  # result of calls that resulted in decrease
        self.list_callback_inp = []  # only appends inputs on callback, as such they correspond to the iterations
        self.list_callback_res = []  # only appends results on callback, as such they correspond to the iterations

    def load_from_optimizer_dict(self):
        """ Load configuration of the optimizer from the dictionary

        Returns
        -------

        """

        for combinationI in self.optimizerDict['parameter_combinations']:
            self.parameter_combinations.append(combinationI)

    def minimize(self, func_to_minimize, args=None, callback=None):
        """ Loop through the parameter combination list

        Each item in the parameter combination list is a dictionary for the Manning's n of each materialID.

        Parameters
        ----------
        func_to_minimize
        parameter_combination_list
        callback

        Returns
        -------

        """

        #get the material ID list and material name list
        materialID_list, materialName_list = args[0], args[1]

        #total number of combinations
        nCombinations = len(self.parameter_combinations)

        #loop over all combinations in the list
        for idx, combinationI in enumerate(self.parameter_combinations):
            #assemble the ManningNs
            ManningNs = [None] * len(materialID_list)

            #loop over all material IDs in the current combination
            for matID, ManningN in combinationI.items():
                if int(matID) in materialID_list:
                    ManningNs[materialID_list.index(int(matID))] = ManningN
                else:  #the current material ID is not in the list; something is wrong
                    raise Exception("The material ID in the enumerator does not match with the list of calibration parameters.")

            #also need to check whether all materials have been dealt with
            for ManningN in ManningNs:
                if ManningN is None:
                    raise Exception("Not all materials are specified in the enumertor.")


            #print information to screen
            print("Parameter combination #", idx+1, "/ ", nCombinations)
            print("     Material names: ", materialName_list)
            print("     Material IDs:   ", materialID_list)
            print("     Manning Ns:     ", ManningNs)

            #call the function to be minimized
            score = func_to_minimize(ManningNs, materialID_list, materialName_list)

            #record the best combination and the corresponding calibration score
            if idx == 0: #first combination
                self.best_combination_index = idx
                self.best_calibration_score = score
            elif self.best_calibration_score > score:
                self.best_combination_index = idx
                self.best_calibration_score = score

            #call the callback
            callback(ManningNs)

        #report the best calibration result
        print("The best calibration parameter combination is # ", self.best_combination_index+1,  #report in 1-based
              "with a minimum calibration score of ", self.best_calibration_score)
        print("The best parameter values are: ")
        print("     Material names: ", materialName_list)
        print("     Material IDs:   ", materialID_list)
        print("     Manning Ns:     ", self.parameter_combinations[self.best_combination_index])  #python is 0-based


    def write_optimization_results_to_csv(self, parameterNames=[]):
        """ Write optimization results to files in csv format

        Returns
        -------

        """

        #write out the calibration parameters' values and calibration errors for all iterations
        fid = open(self.type+"_"+"calibration_results"+".csv", "w")

        #loop through all parameters
        for parameterI in range(len(self.list_callback_inp[0])):
            if len(parameterNames) == 0:
                fid.write("Parameter-%d," % parameterI)
            else:
                fid.write("%s," % parameterNames[parameterI])

        fid.write("Calibration-Error\n")

        #loop through all iterations
        for iterationI in range(len(self.list_callback_inp)):
            for parameterI in range(len(self.list_callback_inp[0])):
                fid.write("%10.5e," % self.list_callback_inp[iterationI][parameterI])

            fid.write("%10.5e\n" % self.list_callback_res[iterationI])

        fid.close()