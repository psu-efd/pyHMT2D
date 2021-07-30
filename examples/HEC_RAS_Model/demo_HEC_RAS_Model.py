"""
Test the HEC_RAS_Model class

It demonstrates how to use pyHMT2D to control the run of HEC-RAS model.
"""

import pyHMT2D

def run_HEC_RAS():

    #create a HEC-RAS model instance
    #my_hec_ras_model = pyHMT2D.RAS_2D.HEC_RAS_Model(version="5.0.7",faceless=False)
    my_hec_ras_model = pyHMT2D.RAS_2D.HEC_RAS_Model(version="6.0.0", faceless=False)

    #initialize the HEC-RAS model
    my_hec_ras_model.init_model()

    print("Hydraulic model name: ", my_hec_ras_model.getName())
    print("Hydraulic model version: ", my_hec_ras_model.getVersion())

    #open a HEC-RAS project
    my_hec_ras_model.open_project("Muncie2D.prj", "Terrain/TerrainMuncie_composite.tif")

    #run the HEC-RAS model's current project
    my_hec_ras_model.run_model()

    #close the HEC-RAS project
    my_hec_ras_model.close_project()

    #quit HEC-RAS
    my_hec_ras_model.exit_model()

def convert_HEC_RAS_to_VTK():
    """ Convert HEC-RAS result to VTK

    Returns
    -------

    """

    my_ras_2d_data = pyHMT2D.RAS_2D.RAS_2D_Data("Muncie2D.p01.hdf",
                                         "Terrain/TerrainMuncie_composite.tif")

    my_ras_2d_data.saveHEC_RAS2D_results_to_VTK(lastTimeStep=True)


if __name__ == "__main__":

    #run_HEC_RAS()

    convert_HEC_RAS_to_VTK()

    print("All done!")

