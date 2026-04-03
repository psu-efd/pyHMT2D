import os
import pyHMT2D
from pyHMT2D.Hydraulic_Models_Data.RAS_2D.HEC_RAS_Model import HEC_RAS_Project

def main():

    #Go into the HEC-RAS directory
    os.chdir("HEC-RAS")

    #build the RAS_to_SRH_Converter object using the new API
    my_ras_to_srh_converter = pyHMT2D.Misc.RAS_to_SRH_Converter("Muncie2D.prj",
                                                                "p01",
                                                                "Muncie")

    #convert to SRH-2D (mesh and material) and save
    print("Convert HEC-RAS 2D mesh and material to SRH-2D.")
    my_ras_to_srh_converter.convert_to_SRH()

    #also convert HEC-RAS 2D result data to VTK (for comparision with SRH-2D)
    print("Convert HEC-RAS 2D results to VTK.")
    project = HEC_RAS_Project("Muncie2D.prj")
    plan = project.get_plan("p01")
    my_ras_2d_data = plan.load_results()

    my_ras_2d_data.saveHEC_RAS2D_results_to_VTK(lastTimeStep=True)

    #Go back to the original directory
    os.chdir("..")


if __name__ == "__main__":
    main()

    print("All done!")
