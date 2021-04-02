"""
A Python class to convert HEC-RAS case, mainly mesh and material (Manning's n), to SRH-2D case
"""

import pyHMT2D

class RAS_to_SRH_Converter:
    """A Python class to convert HEC-RAS case, mainly mesh and material (Manning's n), to SRH-2D case

    Typical work flow is as follows:
    1. create 2D case in HEC-RAS
    2. call RAS to SRH converter: create srhgeom and srhmat files
    3. manually create and adjust srhhydro file
    4. run SRH-2D

    Attributes
    ---------
    RASPlanResultFileName : str
        HEC-RAS result from a plan, such as caseName.p01.hdf
    RASTerrainGeoTiffFileName : str
        HEC-RAS Terrain's GeoTiff file name
    SRHCaseName : str
        Case name for SRH-2D; the resulted files will be SRHCaseName.srhgeom, SRHCaseName.srhmat
    """

    def __init__(self, RASPlanResultFileName, RASTerrainGeoTiffFileName, SRHCaseName):
        self.RASPlanResultFileName = RASPlanResultFileName  # HEC-RAS result from a plan, such as caseName.p01.hdf
        self.RASTerrainGeoTiffFileName = RASTerrainGeoTiffFileName #HEC-RAS Terrain's GeoTiff file name
        self.SRHCaseName = SRHCaseName

    def convert_to_SRH(self):
        """Convert the case to SRH

        Returns
        -------


        """
        #create the RAS_2D_Data object
        ras_2d_data_obj = pyHMT2D.RAS_2D.RAS_2D_Data(self.RASPlanResultFileName, self.RASTerrainGeoTiffFileName)

        ras_2d_data_obj.exportSRHGEOMFile(self.SRHCaseName)

        ras_2d_data_obj.exportSRHMATFile(self.SRHCaseName)