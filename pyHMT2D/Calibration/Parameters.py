
import pyHMT2D

class Parameter(object):
    """ Calibration parameter base class


    Attributes
    ----------

    """

    def __init__(self, active=True):
        """Parameter class constructor

        Parameters
        ----------
        active : bool, optional
            whether this parameter is active or not (default: False). Active refers to whether the parameter
            participate in calibration or scenarios building.


        """

        #type
        self.type = "Parameter_Base"

        #name
        self.name = "Parameter Base"

        #active or not (default: True)
        self.active = True

        self.min = 0.0
        self.max = 0.0

class Parameter_ManningN(Parameter):
    """ ManningN parameter class

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

        self.name = self.type + "_" + str(self.materialID)

        # material name
        self.material_name = parameterDict["material_name"]

        # initial guess for the optimizer
        if "initial_guess" in parameterDict:
            self.initial_guess = parameterDict["initial_guess"]
        else:
            self.initial_guess = None

        # value for the Manning's n
        if "value" in parameterDict:
            self.value = parameterDict["value"]
        else:
            self.value = None

        # bound: min and max
        self.min = parameterDict["min"]
        self.max = parameterDict["max"]

        # some sanity check:
        if self.min > self.max:
            raise Exception("In ManningN parameter specification, the min and larger than the max.")

        if self.initial_guess is not None:
            if self.initial_guess < self.min or \
                    self.initial_guess > self.max:
                raise Exception("In ManningN parameter specification, the initial guess is not in the [min, max] range. Check. Exiting ...")

        if "active" in parameterDict:
            if parameterDict["active"] == "True":
                self.active = True
            elif parameterDict["active"] == "False":
                self.active = False
            else:
                raise  Exception("In Parameter dictionary, \"active\" should be either \"True\" or \"False\". "
                                 "Please check. Exiting ...")

class Parameter_InletQ(Parameter):
    """ InletQ parameter class

    Attributes
    ----------
    bcID : int
        bcID for the InletQ
    value: float
        value for the InletQ
    initial_guess: float
        initial guess of the InletQ
    min, max: float
        min and max of the calibration parameter bound

    """

    def __init__(self, parameterDict):
        """Parameter_InletQ class constructor

        Parameters
        ----------
        parameterDict : dict
            dictionary that contains InletQ parameter information.

        """

        Parameter.__init__(self)

        # type
        self.type = "InletQ"

        # bcID of the InletQ
        self.bcID = parameterDict["bcID"]

        self.name = self.type + "_" + str(self.bcID)

        # bc name
        self.bc_name = parameterDict["bc_name"]

        # initial guess for the optimizer
        if "initial_guess" in parameterDict:
            self.initial_guess = parameterDict["initial_guess"]
        else:
            self.initial_guess = None

        # value for the Manning's n
        if "value" in parameterDict:
            self.value = parameterDict["value"]
        else:
            self.value = None

        # bound: min and max
        self.min = parameterDict["min"]
        self.max = parameterDict["max"]

        # some sanity check:
        if self.min > self.max:
            raise Exception("In InletQ parameter specification, the min and larger than the max.")

        if self.initial_guess is not None:
            if self.initial_guess < self.min or \
                    self.initial_guess > self.max:
                raise Exception("In InletQ parameter specification, the initial guess is not in the [min, max] range. Check. Exiting ...")

        if "active" in parameterDict:
            if parameterDict["active"] == "True":
                self.active = True
            elif parameterDict["active"] == "False":
                self.active = False
            else:
                raise  Exception("In Parameter dictionary, \"active\" should be either \"True\" or \"False\". "
                                 "Please check. Exiting ...")

class Parameter_ExitH(Parameter):
    """ ExitH parameter class

    Attributes
    ----------
    bcID : int
        bcID for the ExitH
    value: float
        value for the ExitH
    initial_guess: float
        initial guess of the ExitH
    min, max: float
        min and max of the calibration parameter bound

    """

    def __init__(self, parameterDict):
        """Parameter_ExitH class constructor

        Parameters
        ----------
        parameterDict : dict
            dictionary that contains ExitH parameter information.

        """

        Parameter.__init__(self)

        # type
        self.type = "ExitH"

        # bcID of the ExitH
        self.bcID = parameterDict["bcID"]

        self.name = self.type + "_" + str(self.bcID)

        # bc name
        self.bc_name = parameterDict["bc_name"]

        # initial guess for the optimizer
        if "initial_guess" in parameterDict:
            self.initial_guess = parameterDict["initial_guess"]
        else:
            self.initial_guess = None

        # value for the Manning's n
        if "value" in parameterDict:
            self.value = parameterDict["value"]
        else:
            self.value = None

        # bound: min and max
        self.min = parameterDict["min"]
        self.max = parameterDict["max"]

        # some sanity check:
        if self.min > self.max:
            raise Exception("In ExitH parameter specification, the min and larger than the max.")

        if self.initial_guess is not None:
            if self.initial_guess < self.min or \
                    self.initial_guess > self.max:
                raise Exception("In ExitH parameter specification, the initial guess is not in the [min, max] range. Check. Exiting ...")

        if "active" in parameterDict:
            if parameterDict["active"] == "True":
                self.active = True
            elif parameterDict["active"] == "False":
                self.active = False
            else:
                raise  Exception("In Parameter dictionary, \"active\" should be either \"True\" or \"False\". "
                                 "Please check. Exiting ...")


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

        # sample in the parameter space
        self.sample()


    def build_parameter_list(self):
        """ Build parameter_list

        Returns
        -------

        """

        # loop through every calibration parameter in the dictionary
        for parameterDict in self.parametersDict:
            if parameterDict["type"] == "ManningN":
                #construct the Parameter object
                currParameter = Parameter_ManningN(parameterDict)

                #append the Parameter_ManningN object to the list
                self.parameter_list.append(currParameter)
            elif parameterDict["type"] == "InletQ":
                #construct the Parameter object
                currParameter = Parameter_InletQ(parameterDict)

                #append the Parameter_InletQ object to the list
                self.parameter_list.append(currParameter)
            elif parameterDict["type"] == "ExitH":
                #construct the Parameter object
                currParameter = Parameter_ExitH(parameterDict)

                #append the Parameter_ExitH object to the list
                self.parameter_list.append(currParameter)
            else:
                raise Exception("The specified calibration parameter type", parameterDict["type"],
                                "is currently not supported. Support type is ManningN. Exiting ...")

    def get_parameter_list(self):
        if len(self.parameter_list)==0:
            self.build_parameter_list()
            return self.parameter_list

        return self.parameter_list

    def get_parameter_space(self):
        """
        Build and return the parameter space.
        [
          (min, max),   #min and max for each parameter in the oder of the configuration list
          (min, max),
          ...
        ]

        Returns
        -------

        """

        parameter_space_list = []

        for parameter in self.parameter_list:
            parameter_space_list.append((parameter.min, parameter.max))

        return parameter_space_list

    def get_parameter_name_list(self):
        """
        Build and return the list of parameter names

        Returns
        -------

        """

        parameter_name_list = []

        for parameter in self.parameter_list:
            parameter_name_list.append(parameter.name)

        return parameter_name_list


    def get_ManningN_Info_list(self):
        """ Get the ManningN's information

        This function returns the following:
        1. a list of calibration materialIDs
        2. a list of corresponding material names
        3. a list of corresponding Manning's n initial guess values
        4. a list of corresponding Manning's n minimum value
        5. a list of corresponding Manning's n maximum value

        Returns
        -------
        materialID_list : list
            a list of calibration materialIDs
        materialName_list : list
            a list of calibration material name
        initial_guess_list : list
            a list of Manning's n initial guesses
        ManningN_min_list : list
            a list of Manning's n minimum values
        ManningN_max_list : list
            a list of Manning's n maximum values


        """

        materialID_list = []
        materialName_list = []
        initial_guess_list = []
        ManningN_min_list = []
        ManningN_max_list =[]

        #loop through all calibration parameter objects (whose type is "ManningN")
        for parameterI in self.parameter_list:
            if parameterI.type == "ManningN" and parameterI.active:
                materialID_list.append(parameterI.materialID)
                materialName_list.append(parameterI.material_name)
                initial_guess_list.append(parameterI.initial_guess)
                ManningN_min_list.append(parameterI.min)
                ManningN_max_list.append(parameterI.max)

        #check the uniqueness of elements in materialID_list. If there are
        #repeated materialIDs, something is not right.
        if len(materialID_list) > len(set(materialID_list)):
            raise Exception("Manning's n calibration materailIDs are not unique. Please check. Exiting ...")

        return materialID_list, materialName_list, initial_guess_list, ManningN_min_list, ManningN_max_list

    def sample(self):
        """
        Sample in the parameter space.



        Returns
        -------

        """
