class HydraulicData(object):
    """ Hydraulic data base class

    Attributes
    ----------
        name : str
            name of the hydraulic mode that this data is for, e.g., SRH-2D, HEC-RAS,
            or just the name of the data set, e.g., terrain, flow

    """

    def __init__(self, name):
        self.name = name

    def getName(self):
        return self.name

    def __str__(self):
        return "Data name: %s " % (self.name)

    def init_data(self):
        """ Initialize the data

        Returns
        -------

        """
        pass

    def reset_data(self):
        """ Reset the data

        Returns
        -------

        """
        pass

    def clear_data(self):
        """ clear the data

        Returns
        -------

        """
        pass

    def modify_ManningsN(self, materialIDs, newManningsNValues, materialNames):
        """Modify materialIDs's Manning's n values to new values

        Parameters
        ----------
        materialID : int
            material ID
        newManningsNValue : float
            new Manning's n value
        materialName : str
            name of the material

        Returns
        -------

        """

        pass