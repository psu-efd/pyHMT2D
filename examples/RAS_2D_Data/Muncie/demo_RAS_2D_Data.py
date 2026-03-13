import pyHMT2D
import shutil
from pathlib import Path

def demo_RAS_2D_Data():
    """ Demonstrate the use of the RAS_2D_Data class to process RAS 2D data

    Returns
    -------

    """

    print("Demonstrating the use of the RAS_2D_Data class to process RAS 2D data...")

    #create a RAS_2D_Data object
    my_ras_2d_data = pyHMT2D.RAS_2D.RAS_2D_Data("Muncie2D.p01.hdf", "Terrain/TerrainMuncie_composite.tif")
    #my_ras_2d_data = pyHMT2D.RAS_2D.RAS_2D_Data("Muncie2D.p01.hdf")  #not working because this case uses a composite terrain file (not sure how to handle that yet)

    #save RAS 2D result to VTK (only last time step in this example; see more options in the documentation for
    #the saveHEC_RAS2D_results_to_VTK(...) function.
    my_ras_2d_data.saveHEC_RAS2D_results_to_VTK(lastTimeStep=True)


def copy_muncie_case_files():
    """Copy all Muncie2D case files/directories into the current directory."""

    print("Copying Muncie2D case files to the current directory...")
    src = Path("../../assets/RAS_models/Muncie2D")
    dst = Path(".")

    for item in src.iterdir():
        target = dst / item.name
        if item.is_dir():
            shutil.copytree(item, target, dirs_exist_ok=True)
        else:
            shutil.copy2(item, target)

if __name__ == "__main__":

    # copy the Muncie2D case files to the current directory
    copy_muncie_case_files()

    demo_RAS_2D_Data()

    print("All done!")