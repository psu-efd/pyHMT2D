
from ..Misc import json_dict_type_correction

import numpy as np

np.random.seed(123)

from skopt.space import Space
from skopt.sampler import Sobol
from skopt.sampler import Lhs
from skopt.sampler import Halton
from skopt.sampler import Hammersly
from skopt.sampler import Grid

class Sampler(object):
    """ Sampler base class


    Attributes
    ----------

    """

    def __init__(self, samplerDict):
        """Optimizer class constructor

        Parameters
        ----------
        samplerDict : dict
            dictionary that contains information on the sampler.

        """

        self.type = "Sampler_Base"

        self.samplerDict = samplerDict

    def load_from_sampler_dict(self):
        """ Load configuration of the sampler from dictionary

        Returns
        -------

        """

        pass

    def generate(self, space, n_samples):
        """
        Generate and return samples

        Parameters
        ----------
        space: Space object
        n_samples: int
            number of samples to draw

        Returns
        -------

        """

        pass


class Sampler_SkoptSampler(Sampler):
    """ Sampler using skopt.sampler class

    Attributes
    ----------

    """

    def __init__(self, samplerDict):
        """Sampler_SkoptSampler class constructor

        Parameters
        ----------

        """

        Sampler.__init__(self, samplerDict)

        # type
        self.type = "SkoptSampler"

        # method: LHS, Sobol, Halton, and Hammersly
        self.method = None

        # dictionary for the options of the particular method
        self.options = {}

        # load sampler configuration from dictionary
        self.load_from_sampler_dict()

        # create the sampler
        self.lhs = None

        self.create_sampler()

    def load_from_sampler_dict(self):
        """ Load configuration of the sampler from the dictionary

        Returns
        -------

        """

        #get the method: e.g., lhs
        self.method = self.samplerDict['method']

        #load options
        self.lhs_type = "classic"
        if "lhs_type" in self.samplerDict:
            self.lhs_type = self.samplerDict["lhs_type"]

        self.criterion = "maximin"
        if "criterion" in self.samplerDict:
            self.criterion = self.samplerDict["criterion"]

        self.iterations = 1000
        if "iterations" in self.samplerDict:
            self.iterations = self.samplerDict["iterations"]

    def create_sampler(self):
        """
        Create the sampler

        Returns
        -------

        """

        if self.method == "lhs":
            self.lhs = Lhs(lhs_type=self.lhs_type, criterion=self.criterion, iterations=self.iterations)
        else:
            raise Exception("The specified sampling method ", self.method, "is not supported.")

    def generate(self, space, n_samples):
        """
        Generate and return samples

        Parameters
        ----------
        space: Space object
        n_samples: int
            number of samples to draw

        Returns
        -------

        array[n_samples, n_parameters]

        """

        return self.lhs.generate(space.dimensions, n_samples)



class Sampler_Enumerator(Sampler):
    """ Sampler using user-provided enumeration of parameter combinations class

    User needs to provide a list of parameter combinations.

    Attributes
    ----------

    """

    def __init__(self, samplerDict):
        """Sampler_Enumerator class constructor

        Parameters
        ----------

        """

        Sampler.__init__(self, samplerDict)

        # type
        self.type = "enumerator"

        # method: default (not relevant).
        self.method = "None"

        #options:

        # dictionary for the options of the particular method
        self.options = {}

        #parameter_combinations (a list of lists)
        #For each list, it contains values for the paramters in the ["parametric_study"]["parameters"] (order is respected)
        self.parameter_combinations = []

        # load sampler configuration from dictionary
        self.load_from_sampler_dict()

    def load_from_sampler_dict(self):
        """ Load configuration of the sampler from the dictionary

        Returns
        -------

        """

        for combinationI in self.samplerDict['enumerator']:
            self.parameter_combinations.append(combinationI)

