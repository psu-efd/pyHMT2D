"""
SRH_2D_Model:
"""

from pyHMT2D.__common__ import HydraulicModel
from .SRH_2D_Data import *

from .helpers import *
import sys
import subprocess
import os

import threading
import queue
import time

from alive_progress import alive_bar

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
        self._srh_2d_data = SRH_2D_Data(srhhydro_filename, srhgeom_filename, srhmat_filename)

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

        print("Running SRH-2D case with input file:", case_srhhydro_file_name)

        #p = subprocess.run([cmd, case_srhhydro_file_name], stdout=subprocess.PIPE,
        #                   stderr=subprocess.PIPE)

        #p = subprocess.run([cmd, case_srhhydro_file_name], stdout=subprocess.PIPE,
        #                   stderr=subprocess.PIPE, start_new_session=True)

        #print(p.stdout)
        #print(p.stderr)

        #if str.encode("successfully executed") in p.stdout:
        #    print("SRH-2D was successfully done!")
        #else:
        #    print("SRH-2D was not successfully done! Check the output files.")

        #this can get real-time output from the external program like srh_2d.exe
        #with subprocess.Popen([cmd, case_srhhydro_file_name], shell=True, stdout=subprocess.PIPE,
        #                   stderr=subprocess.PIPE, bufsize=1) as sp:
        #    for line in sp.stdout:
        #        print(line)

        #r = Runner([cmd, case_srhhydro_file_name])

        #for stdout, stderr in r.start():
        #    print("STDOUT", stdout)
        #    print("STDERR", stderr)

        try:
            print("Remove existing case_INF.DAT")
            os.remove('backwater_curve_INF.DAT')
        except:
            print("Error while deleting file backwater_curve_INF.DAT")
            sys.exit()

        p = subprocess.Popen([cmd, case_srhhydro_file_name], shell=True, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, universal_newlines=True)

        #while True:
            #output = p.stdout.readline()  #this is blocking
        #    if p.poll() is not None:
        #        break
            #if output:
            #    print(output.strip())

        #    time.sleep(2)

        #    info_data = np.genfromtxt("backwater_curve_INF.DAT",skip_header=1)
        #    print(info_data[-1,0])

        with alive_bar(8600) as bar:  # default setting
            for i in range(100):
                bar()  # call after consuming one item


    def run_pre_model(self):
        """Run SRH-2D Pre to preprocess the simulation case

        Returns
        -------

        """

        cmd = self._srh_pre_path
        case_srhhydro_file_name = self._srh_2d_data.srhhydro_filename

        print("Running SRH-2D Pre with case SRHHYDRO input file:", case_srhhydro_file_name)

        p = subprocess.run([cmd, '3', case_srhhydro_file_name], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)

        #print(p.stdout)
        #print(p.stderr)

        if str.encode("successfully executed") in p.stdout:
            print("SRH-2D Pre was successfully done!")
        else:
            print("SRH-2D Pre was not successfully done! Check the output files.")

    def exit_model(self):
        """Exit the model (SRH-2D specific)

        Returns
        -------

        """

        pass


class Runner(object):

    def __init__(self, cmd: []):
        self.cmd = cmd
        self.return_code = None
        self.process = None # type: subprocess.Popen
        self.run_time = 0


    @property
    def _default_popen_kwargs(self):
        return {
            "env": os.environ.copy(),
            "stdout": subprocess.PIPE,
            "stderr": subprocess.PIPE,
            "shell": True,
            "universal_newlines": True,
            "bufsize": 1,
        }


    def _watch_output(self, process: subprocess.Popen, queue):
        for line in iter(process.stderr.readline, ""):
            queue.put(line)
            if process.poll() is not None:
                return

    @property
    def stdout(self):
        return self.process.stdout

    @property
    def stderr(self):
        return self.process.stderr


    def start(self, wait_limit = 15):

        start_time = time.time()

        pargs = self._default_popen_kwargs
        #if self.cmd is not None:
        #    pargs['cmd'] = self.cmd

        self.process = subprocess.Popen(self.cmd,**pargs)
        self.returned = None
        last_output = time.time()
        q = queue.Queue()

        t = threading.Thread(target=self._watch_output, args=(self.process, q,))
        t.daemon = True
        t.start()

        while self.returned is None:
            self.returned = self.process.poll()

            delay = last_output-time.time()
            if self.returned is None:
                stdout = f"{last_output-time.time()} waited"
                try:
                    stderr = q.get_nowait()
                except queue.Empty:
                    time.sleep(1)
                else:
                    yield stdout, stderr
                    last_output = time.time()

            if delay > wait_limit:
                print("Waited 15 seconds, breaking")
                break

        self.run_time = time.time() - start_time