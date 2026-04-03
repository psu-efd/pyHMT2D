"""
Test the HEC_RAS_Model class

It demonstrates how to use pyHMT2D to control the run of HEC-RAS model.
"""

#if run in the cloud, need to specify the location of pyHMT2D. If pyHMT2D is installed
#with regular pip install, then the following is not necessary.
#import sys
#sys.path.append(r"C:\Users\Administrator\python_packages\pyHMT2D")

import sys
#sys.path.append((r"C:\Users\xzl123\Dropbox\PycharmProjects\pyHMT2D"))

#print('\n'.join(sys.path))

import pyHMT2D
import shutil
from pathlib import Path
from pyHMT2D.Hydraulic_Models_Data.RAS_2D.HEC_RAS_Model import HEC_RAS_Project

def run_HEC_RAS():

    print("Running HEC-RAS model ...")

    #create a HEC-RAS model instance
    my_hec_ras_model = pyHMT2D.RAS_2D.HEC_RAS_Model(version="6.6", faceless=False)

    #initialize the HEC-RAS model
    my_hec_ras_model.init_model()

    print("Hydraulic model name: ", my_hec_ras_model.getName())
    print("Hydraulic model version: ", my_hec_ras_model.getVersion())

    #open a HEC-RAS project
    my_hec_ras_model.open_project("Muncie2D.prj")

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

    print("Converting HEC-RAS result to VTK...")

    project = HEC_RAS_Project("Muncie2D.prj")
    plan = project.get_plan("p01")
    my_ras_2d_data = plan.load_results()

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

    run_HEC_RAS()

    convert_HEC_RAS_to_VTK()

    print("All done!")
