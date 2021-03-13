import numpy as np


class HydraulicModel(object):
    """ Hydraulic model base class

    Attributes:
        name: {string} -- name of the hydraulic model, e.g., SRH-2D, HEC-RAS
        version: {string} -- version of the hydaulic model, e.g., 3.3, 5.0.7, 6.0
        init_model: {function} -- initialize the model
        run_model: {function} -- run the model
        end_model: {function} -- end the model
    """
    def __init__(self, name, version):
        self.name = name
        self.version = version

    def getName(self):
        return self.name

    def getVersion(self):
        return self.version

    def __str__(self):
        return "%s version %s" % (self.name, self.version)

    def init_model(self):
        """ Initialize the model

        Returns
        -------

        """
        pass

    def run_model(self):
        """ Run the model

        Returns
        -------

        """
        pass

    def exit_model(self):
        """ Exit the model

        Returns
        -------

        """
        pass

def yes_or_no(question):
    """ a utility tool to ask the user yes/no question

    Parameters
    ----------
    question: the text for the question

    Returns
    -------
    True or False

    """
    while True:
        answer = input(question + ' (y/n): ').lower().strip()
        if answer in ('y', 'yes', 'n', 'no'):
            return answer in ('y', 'yes')
        else:
            print('You must answer y, yes, n, or no.')
