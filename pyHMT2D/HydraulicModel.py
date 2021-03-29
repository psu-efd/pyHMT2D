class HydraulicModel(object):
    """ Hydraulic model base class

    Attributes
    ----------
        name : str
            name of the hydraulic model, e.g., SRH-2D, HEC-RAS
        version : str
            version of the hydaulic model, e.g., 3.3, 5.0.7, 6.0
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