"""
SRH_2D_Model:
"""

from pyHMT2D.__common__ import HydraulicModel
from pyHMT2D.Misc.tools import printProgressBar

from .SRH_2D_Data import *

from .helpers import *
import sys
import subprocess
import os
import time

#from alive_progress import alive_bar

class SRH_2D_Model(HydraulicModel):
    """SRH_2D Model

    SRH_2D_Model controls the run of SRH_2D. A case can be loaded and executed.

    Currently SRH_2D_Model has very limited capability to modify geometry, mesh, and flow data.

    Attributes:
        _faceless: {bool} -- whether show the SRH-2D window when it is running
        _srh_path: path to the SRH_2D program (e.g., srh_2d.exe depending on the version)
        _srh_pre_path: path to the SRH_2D preprocessing program (e.g., srh_2d_pre.exe depending on the version)

        _project_file_name: HEC-RAS project name (including the full path)
        _project: HEC_RAS_Project object corresponding to current project

    """
    def __init__(self, version, srh_pre_path, srh_path, extra_dll_path, faceless=False):
        HydraulicModel.__init__(self, "SRH-2D", version)

        #whether run HEC-RAS without GUI interface showing up
        self._faceless = faceless

        #path for srh_2d_pre.exe and srh_2d.exe (names could vary)
        self._srh_pre_path = srh_pre_path
        self._srh_path = srh_path

        #path for library Dlls, such as XMDL, zlib, hdf5, and msvcr used by SRH-2D.
        self._extra_dll_path = extra_dll_path

        #SRH_2D_Data object to hold information about simulation case
        #including srhhydro, srhgeom, srhmat, etc.
        self._srh_2d_data = None

    def init_model(self):
        """Initialize SRH-2D model

        Returns
        -------

        """
        print("Initializing SRH-2D ...")

        #check whether the specified srh_2d_pre.exe and srh_2d.exe exist
        if not path.isfile(self._srh_pre_path):
            print("The SRH-2D Pre executable file", self._srh_pre_path, "does not exists. Exiting ...")
            sys.exit()

        if not path.isfile(self._srh_path):
            print("The SRH-2D executable file", self._srh_path, "does not exists. Exiting ...")
            sys.exit()

        #check whether the specified extra Dll path exists
        if not path.isdir(self._extra_dll_path):
            print("The extra DLL library path", self._extra_dll_path, "does not exists. Exiting ...")
            sys.exit()

        #add the extra Dll path to "PATH" environment (only temporarily)
        extra_path = os.path.abspath(self._extra_dll_path)
        os.environ['PATH'] = extra_path + ";" + os.environ['PATH']

    def open_project(self, srhhydro_filename, srhgeom_filename, srhmat_filename):
        """Open a SRH-2D project with the specified srhhydro, srhgeom, and srhmat files

        A SRH_2D_Data object will be created to hold the project information

        Parameters
        ----------
        srhhydro_filename: srhhydro file name
        srhgeom_filename: srhgeom file name
        srhmat_filename: srhmat file name

        Returns
        -------

        """
        self._srh_2d_data = SRH_2D_Data(srhhydro_filename)

    def set_simulation_case(self, srh_2d_data):
        """Set the simulation case to srh_2D_data (if it has been created already)

        Parameters
        ----------
        srh_2d_data: an object from class SRH_2D_Data, which should be created before calling

        Returns
        -------

        """
        self._srh_2d_data = srh_2d_data

    def close_project(self):
        """Close the opened SRH-2D project (if any)

        Returns
        -------

        """
        self._srh_2d_data = None

    def run_model(self):
        """Run the SRH-2D model

        Returns
        -------

        """
        cmd = self._srh_path
        case_srhhydro_file_name = self._srh_2d_data.get_case_name() + ".DAT"
        case_INF_file_name = self._srh_2d_data.get_case_name() + "_INF.DAT"

        print("Running SRH-2D case with input file:", case_srhhydro_file_name)

        #get the start and end time in hours
        startTime, endTime = self._srh_2d_data.srhhydro_obj.get_simulation_start_end_time()
        deltaT = self._srh_2d_data.srhhydro_obj.get_simulation_time_step_size()

        print("startTime, endTime (hr) = ", startTime, endTime)
        print("Time step size (s) = ", deltaT)

        try:
            if path.isfile(case_INF_file_name):
                print("Removing case's existing _INF.DAT file ...")
                os.remove(case_INF_file_name)
        except:
            print("Error while deleting case's existing _INF.DAT file. Exiting ...")
            sys.exit()

        p = subprocess.Popen([cmd, case_srhhydro_file_name], shell=True, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)

        #print the simulation progress bar
        totalClicks = 100  # this is the granularity of the progress bar (how many times it needs to be updated)

        while True:
            counter = 0 #counter to track how many times (out of 100) the bar() function has been called

            previous_checked_simulation_time = startTime  #this may not work if it is a restart simulation (?)

            time.sleep(5)  #check the simulation progress every 5 s (this value needs to be updated for each case
            # depending on how long the simulation will be).

            #get the latest simulation information from the INF file
            info_data = np.genfromtxt(case_INF_file_name,skip_header=1)
            #print(info_data[-1,1])  #latest "Time(Hours)"

            time_passed_since_last_check = info_data[-1,1] - previous_checked_simulation_time
            clicks_since_last_check = int(time_passed_since_last_check/(endTime-startTime)*totalClicks)

            counter += clicks_since_last_check

            #print("clicks_since_last_check, counter", clicks_since_last_check, counter)

            if clicks_since_last_check > 0:
                printProgressBar(counter, totalClicks, "Simulation in progress")
                # print("test")   #Don't do this. It will create a new progress bar each time.

            previous_checked_simulation_time = info_data[-1,1]

            if counter > totalClicks:
                break

            if p.poll() is not None:
                break

        print("\n")

        #check whehter the run was successful
        bRunSucessful = False

        for line in p.stdout:
            if str.encode("successfully executed") in line:
                bRunSucessful = True

        if bRunSucessful:
            print("SRH-2D simulation was successfully done!")
        else:
            print("SRH-2D simulation was not successful! Check the output files.")


    def run_pre_model(self):
        """Run SRH-2D Pre to preprocess the simulation case

        Returns
        -------

        """

        cmd = self._srh_pre_path
        case_srhhydro_file_name = self._srh_2d_data.srhhydro_filename

        print("Running SRH-2D Pre with case SRHHYDRO input file:", case_srhhydro_file_name)

        if not path.isfile(case_srhhydro_file_name):
            print("The case SRHHYDRO input file does not exist:", case_srhhydro_file_name, "Exiting ...")
            sys.exit()

        p = subprocess.run([cmd, '3', case_srhhydro_file_name], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)

        if str.encode("successfully executed") in p.stdout:
            print("SRH-2D Pre was successfully done!")
        else:
            print("SRH-2D Pre was not successfully! Check the output files.")

    def exit_model(self):
        """Exit the model (SRH-2D specific)

        Returns
        -------

        """

        pass