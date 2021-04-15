import pyHMT2D

def demo_RAS_2D_Data():
    """ Demonstrate the use of the RAS_2D_Data class to process RAS 2D data

    Returns
    -------

    """

    #create a RAS_2D_Data object
    my_ras_2d_data = pyHMT2D.RAS_2D.RAS_2D_Data("Muncie2D.p01.hdf", "Terrain/TerrainMuncie_composite.tif")

    #save RAS 2D result to VTK (only last time step in this example; see more options in the documentation for
    #the saveHEC_RAS2D_results_to_VTK(...) function.
    my_ras_2d_data.saveHEC_RAS2D_results_to_VTK(lastTimeStep=True)


if __name__ == "__main__":

    demo_RAS_2D_Data()

    print("All done!")