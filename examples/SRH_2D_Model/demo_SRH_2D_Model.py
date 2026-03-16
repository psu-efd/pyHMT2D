"""
Demonstrate the use of SRH_2D_Model class to run a SRH-2D simulation

It demonstrates how to use pyHMT2D to control the run of SRH-2D and convert results to VTK.
"""

import platform
import pyHMT2D

def run_SRH_2D():
    """Run the SRH-2D case."""

    print("Running SRH-2D model ...")

    version = "3.7.1"
    srh_pre_path = r"C:\Program Files\SMS 13.4 64-bit\python\Lib\site-packages\srh2d_exe\SRH_Pre_Console.exe"
    srh_path = r"C:\Program Files\SMS 13.4 64-bit\python\Lib\site-packages\srh2d_exe\SRH-2D_Console.exe"
    extra_dll_path = r"C:\Program Files\SMS 13.4 64-bit\python\Lib\site-packages\srh2d_exe"

    my_srh_2d_model = pyHMT2D.SRH_2D.SRH_2D_Model(
        version, srh_pre_path, srh_path, extra_dll_path, faceless=False
    )
    my_srh_2d_model.init_model()

    print("Hydraulic model name: ", my_srh_2d_model.getName())
    print("Hydraulic model version: ", my_srh_2d_model.getVersion())

    my_srh_2d_model.open_project("Muncie.srhhydro")

    my_srh_2d_model.run_pre_model()
    my_srh_2d_model.run_model()

    my_srh_2d_model.close_project()
    my_srh_2d_model.exit_model()


def convert_SRH_2D_to_VTK():
    """Convert SRH-2D XMDF results to VTK."""

    print("Converting SRH-2D result to VTK...")

    my_srh_2d_data = pyHMT2D.SRH_2D.SRH_2D_Data("Muncie.srhhydro")
    bNodal = False
    xmdf_file = my_srh_2d_data.get_case_name() + "_XMDFC.h5"
    my_srh_2d_data.readSRHXMDFFile(xmdf_file, bNodal)
    my_srh_2d_data.outputXMDFDataToVTK(bNodal, lastTimeStep=True, dir="")


if __name__ == "__main__":

    run_SRH_2D()

    convert_SRH_2D_to_VTK()

    print("All done!")
