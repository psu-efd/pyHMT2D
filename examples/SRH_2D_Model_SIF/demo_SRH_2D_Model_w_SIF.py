"""
Test the SRH_2D_Model class

Run SRH-2D simulation with pyHMT2D
"""
import platform
import pyHMT2D

pyHMT2D.gVerbose = True

def run_SRH_2D(srhcontrol_file, system_name):
    """ Run the SRH-2D case

    Returns
    -------

    """   

    if system_name == "Windows":
        #the following should be modified based on your installation of SRH-2D on Windows
        version = "3.6.5"
        srh_pre_path = r"C:\Program Files\SMS 13.3 64-bit\python\Lib\site-packages\srh2d_exe\SRH_Pre_Console.exe"
        srh_path = r"C:\Program Files\SMS 13.3 64-bit\python\Lib\site-packages\srh2d_exe\SRH-2D_Console_v365.exe"
        extra_dll_path = r"C:\Program Files\SMS 13.3 64-bit\python\Lib\site-packages\srh2d_exe"
    elif system_name == "Linux":
        #the following should be modified based on your installation of SRH-2D on Linux
        version = "3.6.2"
        srh_pre_path = r"../../bin/srh2dpre"
        srh_path = r"../../bin/srh2d"
        extra_dll_path = r"../../bin"
    else:
        raise ValueError("Unsupported operating system: " + system_name)

    #create a SRH-2D model instance
    my_srh_2d_model = pyHMT2D.SRH_2D.SRH_2D_Model(version, srh_pre_path,
                       srh_path, extra_dll_path, faceless=False)

    #initialize the SRH-2D model
    my_srh_2d_model.init_model()

    print("Hydraulic model name: ", my_srh_2d_model.getName())
    print("Hydraulic model version: ", my_srh_2d_model.getVersion())

    #open a SRH-2D project
    my_srh_2d_model.open_project(srhcontrol_file)

    #run SRH-2D Pre to preprocess the case
    my_srh_2d_model.run_pre_model()

    #run the SRH-2D model's current project
    my_srh_2d_model.run_model()

    #close the SRH-2D project
    my_srh_2d_model.close_project()

    #quit SRH-2D
    my_srh_2d_model.exit_model()


def convert_SRH_2D_to_VTK(srhcontrol_file, system_name):
    """ Convert SRH-2D results to VTK


    Returns
    -------

    """

    my_srh_2d_data = pyHMT2D.SRH_2D.SRH_2D_Data(srhcontrol_file)

    #read SRH-2D result in XMDF format (*.h5)
    #wether the XMDF result is nodal or cell center
    bNodal = False

    #read the XMDF file
    #my_srh_2d_data.readSRHXMDFFile("Muncie_XMDFC.h5", bNodal)

    #export the XMDF data to VTK
    #my_srh_2d_data.outputXMDFDataToVTK(bNodal, lastTimeStep=True, dir='')

    #read the SRHC files
    my_srh_2d_data.readSRHCFiles(my_srh_2d_data.srhsif_obj.srhsif_content["Case"])

    #export the SRHC data to VTK
    my_srh_2d_data.outputSRHCDataToVTK(lastTimeStep=True, dir='')

if __name__ == "__main__":


    #determine the operating system
    system_name = platform.system()    

    #test the reading of the SIF file
    #my_sif = pyHMT2D.SRH_2D.SRH_2D_SIF("Muncie_SIF.dat")
    #print(my_sif.srhsif_content)

    #change the Manning's n values
    #my_sif.modify_manning_n(2, 0.012345)

    #save to a new SIF file
    #my_sif.save("Muncie_SIF_new.dat")

    #my_srh_2d_data = pyHMT2D.SRH_2D.SRH_2D_Data("Muncie_SIF.dat")
    #print(my_srh_2d_data.srhsif_obj.srhsif_content["Case"])
    #print(my_srh_2d_data.srhsif_obj.srhsif_content)

    #my_srh_2d_data = pyHMT2D.SRH_2D.SRH_2D_Data("Muncie.srhhydro")
    #print(my_srh_2d_data.srhhydro_obj.srhhydro_content["Case"])
    #print(my_srh_2d_data.srhhydro_obj.srhhydro_content)

    convert_SRH_2D_to_VTK("Muncie_SIF.dat", system_name)


    if system_name == "Windows":        #using srhhydro (on Windows)
        srhcontrol_file = "Muncie.srhhydro"
        #run_SRH_2D(srhcontrol_file, system_name)
        #convert_SRH_2D_to_VTK(srhcontrol_file, system_name)
    elif system_name == "Linux":        #using SIF (on Linux)
        srhcontrol_file = "Muncie_SIF.dat"
        #run_SRH_2D(srhcontrol_file, system_name)
        #convert_SRH_2D_to_VTK(srhcontrol_file, system_name)
    else:
        raise ValueError("Unsupported operating system: " + system_name)

    print("All done!")

