"""
SRH_2D_Model:
"""

from pyHMT2D.Hydraulic_Models_Data import HydraulicModel
from pyHMT2D.Misc.tools import printProgressBar

from .SRH_2D_Data import *

import sys
import subprocess
import os
import time
import platform


class SRH_2D_Model(HydraulicModel):
    """SRH-2D Model which can control simulation runs.

    SRH_2D_Model controls the run of SRH-2D. A case can be loaded and executed. Currently
    SRH_2D_Model has very limited capability to modify geometry, mesh, and flow data. One
    of the main use scenarios is for calibration or Monte Carlo simulations.

    Attributes:
        _faceless : bool
            whether show the SRH-2D window when it is running
        _srh_path : str
            path to the SRH_2D program (e.g., srh_2d.exe depending on the version of SRH-2D and the operating system)
        _srh_pre_path : str
            path to the SRH_2D preprocessing program (e.g., srh_2d_pre.exe depending on the version of SRH-2D and the operating system)
        _extra_dll_path : str
            path for library Dlls, such as XMDL, zlib, hdf5, and msvcr used by SRH-2D (only relevant for Windows).
        _srh_2d_data : SRH_2D_Data
            a SRH_2D_Data object to hold the simulation case data.

    """
    def __init__(self, version, srh_pre_path, srh_path, extra_dll_path, faceless=False):
        """SRH_2D_Model constructor

        Parameters
        ----------
        version : str
            version number of SRH-2D
        srh_pre_path : str
            path to the SRH_2D preprocessing program (e.g., srh_2d_pre.exe depending on the version)
        srh_path : str
            path to the SRH_2D program (e.g., srh_2d.exe depending on the version)
        extra_dll_path : str
            path for library Dlls, such as XMDL, zlib, hdf5, and msvcr used by SRH-2D
        faceless : bool, optional
            whether show the SRH-2D window when it is running
        """

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

        #check the operating system
        self._os_type = platform.system()


    def init_model(self):
        """Initialize SRH-2D model

        It will check the existance of SRH-2D executable files and add extra DLL path
        to the system PATH environment (only relevant for Windows).

        Returns
        -------

        """

        print("Initializing SRH-2D ...")

        #check whether the specified executable files exist
        if not path.isfile(self._srh_pre_path):
            print("The SRH-2D Pre executable file", self._srh_pre_path, "does not exists. Exiting ...")
            sys.exit()


        if not path.isfile(self._srh_path):
            print("The SRH-2D executable file", self._srh_path, "does not exists. Exiting ...")
            sys.exit()

        #check whether the specified extra Dll path exists (only relevant for Windows)
        if self._os_type == "Windows":
            if not path.isdir(self._extra_dll_path):
                print("The extra DLL library path", self._extra_dll_path, "does not exists. Exiting ...")
                sys.exit()

            #add the extra Dll path to "PATH" environment (only temporarily)
            extra_path = os.path.abspath(self._extra_dll_path)
            os.environ['PATH'] = extra_path + ";" + os.environ['PATH']

    def open_project(self, srhcontrol_filename):
        """Open a SRH-2D project with the specified srhcontrol file which can be a srhhydro or SIF file

        A SRH_2D_Data object will be created to hold the project information.

        Parameters
        ----------
        srhcontrol_filename : str
            name of the srhcontrol file which can be a srhhydro or SIF file

        Returns
        -------

        """

        #check whether it is a srhhydro or SIF file
        if not srhcontrol_filename.endswith(".srhhydro") and \
           not ("sif" in srhcontrol_filename.lower()):  #contains "sif" or "SIF"
            print("The SRH-2D control file must be either srhhydro or SIF file. Exiting ...")
            sys.exit()

        #check if the srhcontrol file exists
        if not path.isfile(srhcontrol_filename):

            print("The SRH-2D control file", srhcontrol_filename, "does not exists. Exiting ...")
            sys.exit()


        self._srh_2d_data = SRH_2D_Data(srhcontrol_filename)

        #currently, pyHMT2D only supports srhhydro file on Windows and SIF file on Linux
        if self._os_type == "Windows":
            if self._srh_2d_data.control_type == "SRHHydro":
                pass
            elif self._srh_2d_data.control_type == "SIF":
                print("Currently, pyHMT2D only supports SIF file on Linux, not Windows. Exiting ...")
                sys.exit()
        elif self._os_type == "Linux":
            if self._srh_2d_data.control_type == "SRHHydro":
                print("Currently, pyHMT2D only supports SRHHydro file on Windows, not Linux. Exiting ...")
                sys.exit()
            elif self._srh_2d_data.control_type == "SIF":
                pass


    def set_simulation_case(self, srh_2d_data):
        """Set the simulation case to srh_2D_data (if it has been created already)

        Parameters
        ----------
        srh_2d_data : SRH_2D_Data
            an object from class SRH_2D_Data, which should be created before calling

        Returns
        -------

        """

        self._srh_2d_data = srh_2d_data

    def get_simulation_case(self):
        """Get the simulation case to srh_2D_data

        Parameters
        ----------

        Returns
        -------
        srh_2d_data : SRH_2D_Data
            an object from class SRH_2D_Data, which should be created before calling

        """

        return  self._srh_2d_data

    def close_project(self):
        """Close the opened SRH-2D project (if any)

        It will set _srh_2d_data to None.

        Returns
        -------

        """

        self._srh_2d_data = None

    def run_model(self, sleepTime = 10.0, bShowProgress=True):
        """Run the SRH-2D model

        It will run the current SRH-2D project (case) and show a progress bar. The case
        has to be created before with open_project(...) or set_simulation_case().

        :param sleepTime: float
            sleep time between each check on the progress
        :param bShowProgress: bool
            whether show the progress bar

        Returns
        -------

        """

        cmd = self._srh_path
        case_file_name = self._srh_2d_data.get_case_name() + ".DAT"           #e.g., Cimarron.DAT
        case_INF_file_name = self._srh_2d_data.get_case_name() + "_INF.DAT"   #e.g., Cimarron_INF.DAT

        print("Running SRH-2D case with input file:", case_file_name)

        #get the start and end time in hours
        startTime, endTime = 0.0, 0.0
        deltaT = 0.0
        if self._srh_2d_data.control_type == "SRHHydro":
            startTime, endTime = self._srh_2d_data.srhhydro_obj.get_simulation_start_end_time()
            deltaT = self._srh_2d_data.srhhydro_obj.get_simulation_time_step_size()
        elif self._srh_2d_data.control_type == "SIF":
            startTime, endTime = self._srh_2d_data.srhsif_obj.get_simulation_start_end_time()
            deltaT = self._srh_2d_data.srhsif_obj.get_simulation_time_step_size()

        print("startTime, endTime (hr) = ", startTime, endTime)
        print("Time step size (s) = ", deltaT)

        try:
            if path.isfile(case_INF_file_name):
                print("Removing case's existing _INF.DAT file ...")
                os.remove(case_INF_file_name)
        except:
            print("Error while deleting case's existing _INF.DAT file. Exiting ...")
            sys.exit()

        p = subprocess.Popen([cmd, case_file_name], shell=False, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)

        if bShowProgress:
            # print the simulation progress bar
            totalClicks = 100  # this is the granularity of the progress bar (how many times it needs to be updated)

            while True:
                counter = 0  # counter to track how many times (out of 100) the bar() function has been called

                previous_checked_simulation_time = startTime  # this may not work if it is a restart simulation (?)

                time.sleep(10)  # check the simulation progress every 5 s (this value needs to be updated for each case
                # depending on how long the simulation will be).

                # get the latest simulation information from the INF file
                if os.path.isfile(case_INF_file_name):
                    info_data = np.genfromtxt(case_INF_file_name, skip_header=1)
                    # print(info_data[-1,1])  #latest "Time(Hours)"
                else:
                    continue

                # need to consider the scenarios that there is zero line, one line, and multiple lines of INFO output
                current_simulation_time = startTime
                try:
                    if len(info_data.shape) == 0:
                        current_simulation_time = startTime
                    elif len(info_data.shape) == 1:
                        current_simulation_time = info_data[1]
                    elif len(info_data.shape) == 2:
                        current_simulation_time = info_data[-1, 1]

                    time_passed_since_last_check = current_simulation_time - previous_checked_simulation_time
                    clicks_since_last_check = int(time_passed_since_last_check / (endTime - startTime) * totalClicks)

                    counter += clicks_since_last_check

                    # print("clicks_since_last_check, counter", clicks_since_last_check, counter)

                    if clicks_since_last_check > 0:
                        printProgressBar(counter, totalClicks, "Simulation in progress")
                        # print("test")   #Don't do this. It will create a new progress bar each time.

                    previous_checked_simulation_time = current_simulation_time

                except IndexError:  # try to catch the occational index error if the INF write is not compplete.
                    continue

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

        #p.stdout.close()
        #p.kill()

        return bRunSucessful


    def run_pre_model(self):
        """Run SRH-2D Pre to preprocess the simulation case

        Run SRH-2D's preprocessing step for current case.

        Returns
        -------

        """

        cmd = self._srh_pre_path
        case_srhcontrol_file_name = self._srh_2d_data.srhcontrol_filename

        print("Running SRH-2D Pre with case control input file (srhhydro or SIF):", case_srhcontrol_file_name)

        if not path.isfile(case_srhcontrol_file_name):
            print("The case control file does not exist:", case_srhcontrol_file_name, "Exiting ...")
            sys.exit()

        #depending on the control file type, call the corresponding preprocessing program
        if self._srh_2d_data.control_type == "SRHHydro":
            p = subprocess.run([cmd, '3', case_srhcontrol_file_name], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
        elif self._srh_2d_data.control_type == "SIF":
            #get the base case name, e.g., Cimarron_SIF.dat -> Cimarron
            base_case_name = self._srh_2d_data.get_case_name()
            p = subprocess.run([cmd, base_case_name], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)

        if str.encode("successfully executed") in p.stdout:
            print("SRH-2D Pre was successfully done!")
            #p.stdout.close()
            #p.kill()
            return True
        else:
            print("SRH-2D Pre was not successfully! Check the output files.")
            #p.stdout.close()
            #p.kill()
            return False

    def exit_model(self):
        """Exit the model (SRH-2D specific)

        Returns
        -------

        """

        pass