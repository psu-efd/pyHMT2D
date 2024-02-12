"""
RAS_2D_Model:
"""

#import win32com.client as win32  #moved inside code
import os.path

from pyHMT2D.Hydraulic_Models_Data import HydraulicModel

from pyHMT2D.__common__ import gVerbose

from pyHMT2D.Hydraulic_Models_Data.RAS_2D import RAS_2D_Data

from .helpers import *
import sys

class HEC_RAS_Project(object):
    """ HEC-RAS project class (data from .prj file)

    Attributes:
        title : str
            project title
        currentPlanName : str
            name of the current plan
        geom_file_list : list
            list of all geometry files
        flow_file_list : list
            list of all flow files
        plan_file_list : list
            list of all plan files
    """

    def __init__(self, title='', currentPlanName='', currentPlanFile='', geom_file_list=[],
                 flow_file_list=[], plan_file_list=[], plans=[], terrainFileName=''):
        """HEC_RAS_Project class constructor

        Parameters
        ----------
        title : str, optional
            title
        currentPlanName : str, optional
            current plan's name
        currentPlanFile : str, optional
            current plan's file name
        geom_file_list : str, optional
            list of geometry file names
        flow_file_list : str, optional
            list of flow file names
        plan_file_list : str, optional
            list of plan file names
        plans : str
            list of HEC_RAS_Plan objects
        terrainFileName : str
            name of terrain file (not necessary the terrain used by HEC-RAS; it is used for interpolation)
        """

        self.title = title
        self.currentPlanName = currentPlanName
        self.currentPlanFile = currentPlanFile
        self.geom_file_list = geom_file_list
        self.flow_file_list = flow_file_list
        self.plan_file_list = plan_file_list

        #Plans in the opened HEC-RAS project: a list of HEC_RAS_Plan objects
        self.plans = plans

        self.terrainFileName = terrainFileName

    def set_current_plan_file_name(self, currentPlanName, currentPlanFile):
        """

        Parameters
        ----------
        currentPlanName : str, optional
            current plan's name
        currentPlanFile : str, optional
            current plan's file name

        Returns
        -------

        """

        self.currentPlanName = currentPlanName
        self.currentPlanFile = currentPlanFile

    def get_current_plan_file_name(self):
        return self.currentPlanName, self.currentPlanFile

    def get_terrainFileName(self):
        return self.terrainFileName

    def __del__(self):
        """ Destructor

        Clear out lists and objects (may not be necessary; but just to make sure
        we don't cause memory leak)

        Returns
        -------

        """

        self.geom_file_list.clear()
        self.flow_file_list.clear()
        self.plan_file_list.clear()
        self.plans.clear()

    def __str__(self):
        """Define the string representation of the object as HEC-RAS case's title

        Returns
        -------

        """
        return self.title

    def __repr__(self):
        """Define the object as a string representation

        Returns
        -------

        """
        return 'HEC-RAS Project: title = "{}", ' \
               'current plan name = "{}", ' \
               'geometry file list = "{}", ' \
               'flow file list = "{}", ' \
               'plan file list = "{}", ' \
               'plans = "{}"'.format(self.title, self.currentPlanName,
                                         self.geom_file_list, self.flow_file_list,
                                         self.plan_file_list, self.plans)

class HEC_RAS_Plan(object):
    """ HEC-RAS plan (data from .p## file)

    Attributes:
         title : str
            title of the plan
         steady : bool
            whether it is steady or not
         geom_file : str
            geometry file name, e.g., g01, g02, etc.
         flow_file : str
            flow file name, e.g, u01, u02, etc.
         plan_file : str
            plan file name, e.g., p01, p02, etc.

    """
    def __init__(self, title='', steady=False, geom_file='', flow_file='', plan_file=''):
        """HEC_RAS_Plan class constructor

        Parameters
        ----------
        title : str, optional
            title of the plan
        steady : bool, optional
            whether the plan is steady or not (unsteady)
        geom_file : str, optional
            geometry file name
        flow_file : str, optional
            flow file name
        plan_file : str, optional
            plan file name
        """

        self.title = title
        self.steady = steady
        self.geom_file = geom_file
        self.flow_file = flow_file
        self.plan_file = plan_file

    def __str__(self):
        return self.title

    def __repr__(self):
        return 'HEC-RAS plan: title = "{}", ' \
               'steady = "{}", ' \
               'geometry file = "{}", ' \
               'flow file = "{}"'.format(self.title, self.steady,
                                         self.geom_file, self.flow_file)

class HEC_RAS_Geometry(object):
    """ Data for a HEC-RAS geometry (from .g## file)

    Not implemented yet.
    """

    pass

class HEC_RAS_SteadyFlow(object):
    """ Data for a HEC-RAS steady flow (from .f## file)

    Not implemented yet.
    """

    pass

class HEC_RAS_UnsteadyFlow(object):
    """ Data for a HEC-RAS unsteady flow (from .u## file)

    Not implemented yet.
    """

    pass

class HEC_RAS_Model(HydraulicModel):
    """HEC-RAS Model class

    HEC-RAS Model controls the run of HEC-RAS. A plan can be loaded and executed.
    Currently HEC_RAS_Model has very limited capability to modify HEC-RAS project,
    plan, geometry, and flow data. It is mainly used for calibration and Monte
    Carlo simulations.

    Attributes:
        _faceless : bool
            whether show the HEC-RAS GUI
        _ras_path : str
            path to the HEC-RAS program (Ras.exe)
        _RASController: object
            RASController
        _project_file_name : str
            HEC-RAS project name (including the full path)
        _project : HEC_RAS_Project
            HEC_RAS_Project object corresponding to current project

    """

    def __init__(self, version, faceless=False):
        """HEC_RAS_Model constructor

        Parameters
        ----------
        version : str
            HEC-RAS version
        faceless : bool, optional
            whether to run HEC-RAS without GUI
        """

        HydraulicModel.__init__(self, "HEC-RAS", version)

        #whether run HEC-RAS without GUI interface showing up
        self._faceless = faceless

        #path to HEC-RAS's Ras.exe, e.g., C:\Program Files (x86)\HEC\HEC-RAS\5.0.7\Ras.exe
        self._ras_path = None

        #check the supplied version is supported or not
        if not (version in get_supported_hec_ras_versions()):
            print("Specified HEC-RAS version", version, "is not supported in current version of pyHMT2D.")
            print("Supported HEC-RAS versions are", get_supported_hec_ras_versions())
            print("Exitting.")
            sys.exit()

        self._installed_hec_ras_versions = get_installed_hec_ras_versions() # get installed HEC-RAS versions on the
        # computer
        if (not self._installed_hec_ras_versions): #if the available HEC-RAS list is empty
            print("No HEC-RAS installed on this computer. Install HEC-RAS first. Exitting.")
            sys.exit()
        else:
            print("Available versions of HEC-RAS: ", self._installed_hec_ras_versions)

        #check whether the specified HEC-RAS version is available
        if not (version in get_installed_hec_ras_versions()):
            print("Specified HEC-RAS version", version, "is not installed on this computer. Install it first.",
                                                         "Exitting.")
            sys.exit()

        #RASController
        self._RASController = None

        #HEC-RAS project name (including the full path)
        self._project_file_name = ''

        #An HEC_RAS_Project object for the opened HEC-RAS project
        self._project = None

        #RAS_2D_Data object to hold information about simulation case
        #including mesh, boundary, and results
        self._ras_2d_data = None

    def __del__(self):
        """Some clean up before exit."""

        pass

    def init_model(self):
        """ Initialize a HEC-RAS instance

        Returns
        -------

        """

        print("Initializing HEC-RAS ...")

        import platform

        if platform.system().startswith('Windows'):
            try:
                import win32com.client as win32
            except ImportError:
                raise ImportError('Error in importing pywin32 package. Make sure it has been installed properly.')
        else:
            print("Current OS: ", platform.system())
            raise ImportError('The use of pyHMT2D for HEC-RAS Model is only supported on Windows.')

        # Should we kill all currently running HEC-RAS instances? Probably not a good idea.
        #kill_all_hec_ras()

        # If the HEC-RAS model has already been initialized, do nothing.
        if self._RASController is not None:
            # print and confirm the HEC-RAS version
            print('HEC-RAS version ', self._RASController.HECRASVersion(), 'has already been initialized.')

            return

        # Difference between Dispatch and DispatchEx: the latter launches a new instance
        # Two options: not too much of difference for a user
        # late binding: python does not know what the object (HEC-RAS) can do
        # ras = win32.Dispatch("RAS5x.HECRASController")
        # or
        # early binding: python knows what the object (HEC-RAS) can do

        if (self.getVersion() == '5.0.7'):
            #self._RASController = win32.gencache.EnsureDispatch('RAS507.HECRASController')
            #self._RASController = win32.Dispatch("RAS507.HECRASController")
            self._RASController = win32.DispatchEx("RAS507.HECRASController")

        elif (self.getVersion() == '6.0.0'):
            # As of 03/21, "RAS5x.HECRASController" is the prog_id for HEC-RAS version 6 beta 1 and 2
            # This may change (check in future).
            # Updated 07/21: The latest official release of HEC-RAS v6.0.0's prog_id is "RAS60.HECRASController"
            self._RASController = win32.gencache.EnsureDispatch('RAS60.HECRASController')

        elif (self.getVersion() == '6.1.0'):
            # As of 03/21, "RAS5x.HECRASController" is the prog_id for HEC-RAS version 6 beta 1 and 2
            # This may change (check in future).
            # Updated 07/21: The latest official release of HEC-RAS v6.1.0's prog_id is "RAS60.HECRASController"
            self._RASController = win32.gencache.EnsureDispatch('RAS610.HECRASController')
        elif (self.getVersion() == '6.3.1'):
            self._RASController = win32.gencache.EnsureDispatch('RAS631.HECRASController')
        elif (self.getVersion() == '6.4.1'):
            self._RASController = win32.gencache.EnsureDispatch('RAS641.HECRASController')
        else:
            raise Exception("The specified version of HEC-RAS is not currently supported.")

        # show HEC-RAS main window
        if not self._faceless:
            self._RASController.ShowRas()

        # print and confirm the HEC-RAS version
        print('Successfully launched HEC-RAS version ', self._RASController.HECRASVersion())

        #_ = input("Press ENTER to continue:")


    def open_project(self,projectFileName, terrainFileName):
        """ open the specified HEC-RAS project

        Parameters
        ----------
        projectFileName : str
            project file name including the path
        terrainFileName : str
            file name for the terrain (not necessarily the terrain used by HEC-RAS for the project). This
            terrain file (in GeoTiff format) is used to interpolate the elevation at face points (which HEC-RAS
            does not export to result HDF files). This GeoTiff terrain file can be converted from the terrain used
            in the HEC-RAS project. (TODO: convert in this function internally without user input).

        Returns
        -------
            currentPlanFile : str
                name of current plan file
            currentPlanName: str
                name of current plan

        """

        if gVerbose: print("HEC-RAS opens project: ", projectFileName)

        # check whether the HEC-RAS model has been initialized. If not, call init_model()
        if self._RASController is None:
            self.init_model()

        # check whether the specified projectFileName file exists
        if os.path.isfile(projectFileName):
            #get the absolute path of the project file (it seems HEC-RAS does not like the relative path)
            self._project_file_name = os.path.abspath(projectFileName)
        else:
            error = 'Project file "{}" not found.'.format(projectFileName)

        # open the project
        self._RASController.Project_Open(self._project_file_name)

        title = self._RASController.CurrentProjectTitle()


        ###########################################
        # Build data for the HEC_RAS_Project object
        ###########################################
        PlanCount = 0
        PlanNames = None
        IncludeOnlyPlansInBaseDirectory = False

        PlanCount, PlanNames = self.get_plan_names(IncludeOnlyPlansInBaseDirectory)

        print("There are ", PlanCount, "plan(s) in current project.")
        print('Plan names: ', PlanNames)  # returns plan names

        #infer the current plan from current plan file
        #HEC-RAS does not provide a function to get current plan name
        currentPlanFile = self._RASController.CurrentPlanFile()

        currentPlanName = ''

        for i in range(PlanCount):
            if currentPlanFile == self._RASController.Plan_GetFilename(PlanNames[i])[0]:
                currentPlanName = PlanNames[i]
                break

        print("Current plan name is ", currentPlanName)

        #If plans are not empty, build these plans

        print("Building all the plans in the project ...")

        geom_file_list = []
        flow_file_list = []
        plan_file_list = []

        plans = []

        if PlanCount > 0:
            for i in range(PlanCount):
                print("Plan ", i, ",", PlanNames[i], ":")

                #temporarily set the current plan to be PlanNames[i]; will be restored.
                self._RASController.Plan_SetCurrent(PlanNames[i])

                #get the plan file for the current plan
                plan_file = self._RASController.CurrentPlanFile()
                plan_file_list.append(plan_file)

                #get the geometry file for the current plan
                geom_file = self._RASController.CurrentGeomFile()
                geom_file_list.append(geom_file)

                #get the flow file for the current plan
                flow_file_steady = self._RASController.CurrentSteadyFile()
                flow_file_unsteady = self._RASController.CurrentUnSteadyFile()

                #print("flow_file_steady = ", flow_file_steady)
                #print("flow_file_unsteady = ", flow_file_unsteady)

                #check whether both steady and unsteady flow files are empty
                if (not flow_file_steady) and (not flow_file_unsteady):
                    print("Error: no flow file specified in the current plan. Exitting.")
                    sys.exit()

                steady = False

                if (flow_file_unsteady):
                    print("    It is an unsteady flow plan.")
                    flow_file = flow_file_unsteady
                else:
                    print("    It is a steady flow plan.")
                    steady = True
                    flow_file = flow_file_steady

                flow_file_list.append(flow_file)

                plans.append(HEC_RAS_Plan(PlanNames[i], steady, geom_file, flow_file, plan_file))

        #print(plans)

        #set the current plan back to its original value
        self._RASController.Plan_SetCurrent(currentPlanName)

        #create the HEC_RAS_Project object
        if self._project is not None:
            self._project = None

        self._project = HEC_RAS_Project(title, currentPlanName, currentPlanFile, geom_file_list,
                                        flow_file_list, plan_file_list, plans, terrainFileName)

        #dump _project content to screen
        #print(self._project)

        #create RAS_2D_Data object (can't be here because the plan.hdf file exists only after the plan has been run)
        #self._ras_2d_data = RAS_2D_Data(currentPlanFile + ".hdf", terrainFileName)

        print("Finished building all the plans in the project.")

        return currentPlanName, currentPlanFile

    def set_current_plan(self, planName):
        """

        Parameters
        ----------
        planName : str
            name of the plan to be set as current. The name should be like "Peak_flow_plan", not "p01" or "p02"

        Returns
        -------

        """

        #TODO: should check whether planName is in the list plans of current project

        #Plan_SetCurrent(...) returns a tuple (True/False, planName)
        ret = self._RASController.Plan_SetCurrent(planName)

        if ret[0]:
            print("Call to set_current_plan(...) is successful. Current plan is ", planName)
        else:
            print("Call to set_current_plan(...) failed. Current plan not changed. Check the validity of the plan name: ", planName)

    def get_current_planFile(self):
        """
        Return the current plan file (it seems there is no function to get plan name)

        Returns
        -------

        """

        plan_file = self._RASController.CurrentPlanFile()

        return plan_file

    def load_current_plan_results(self):
        """
        Load the current plan results to ras_2d_data
        Returns
        -------

        """

        currentPlanName, currentPlanFile = self._project.get_current_plan_file_name()
        terrainFileName = self._project.get_terrainFileName()

        #check whether the hdf file exists
        if not os.path.isfile(currentPlanFile):
            raise Exception("The result HDF file for current plan does not exist. Make sure to run HEC-RAS before calling this function.")
        elif not os.path.isfile(terrainFileName):
            raise Exception("The specified terrain file does not exist.")
        else:
            self._ras_2d_data = RAS_2D_Data(currentPlanFile + ".hdf", terrainFileName)


    def get_simulation_case(self, bReload=False):
        """Get the simulation case to RAS_2D_data

        Parameters
        ----------

        Returns
        -------
        ras_2d_data : RAS_2D_Data
            an object from class RAS_2D_Data, which should be created before calling
        bRelaod: bool
            whether to reload the results from the HDF file

        """

        if self._ras_2d_data is None or bReload:
            self.load_current_plan_results()

        return self._ras_2d_data

    def get_current_project(self):
        """Get the current opened project

        Returns
        -------

        """

        if (self._project is not None):
            return self._project
        else:
            raise Exception("HEC-RAS has no project opend. Open the project first.")


    def save_project(self):
        """Save the current project

        Returns
        -------

        """

        if (self._RASController is not None)  and (self._project is not None):
            print("Saving project: ", self._RASController.CurrentProjectTitle())
            self._RASController.Project_Save()

    def close_project(self):
        """Close the current project (if any)

        Returns
        -------

        """

        if (self._RASController is not None) and (self._project is not None):
            print("Closing project: ", self._RASController.CurrentProjectTitle())
            self._RASController.Project_Close()

            self._project = None


    def run_model(self):
        """ Run the HEC-RAS model

        A project and a plan have to be defined and selected before the run.

        Returns
        -------

        """

        #check: HEC-RAS controller has been initialized and a project has been opened.
        if self._RASController is None:
            print("HEC-RAS model has not been initialized yet. Call init_model() first.")
            return

        if not (self._RASController.Project_Current()):
            print("No HEC-RAS project has been opened yet. Call open_project(...) first.")

        print("HEC-RAS is computing the current plan ...")

        if self._faceless:
            self._RASController.Compute_HideComputationWindow()
        else:
            self._RASController.Compute_ShowComputationWindow()

        nmsg = None
        msg = None
        res = self._RASController.Compute_CurrentPlan(nmsg, msg)

        #print computing message
        bRunSucessful = True
        if res[0]:
            print("HEC-RAS computed successfully.")
        else:
            print("HEC-RAS computed unsuccessfully. The HEC-RAS Controller's Compute_CurrentPlan() function returned "
                  "False.")
            bRunSucessful = False

        print("The returned messages are:")
        print("res = ", res)
        print("nmsg = ", nmsg)
        print("msg = ", msg)
        for i in range(res[1]):
            print("    ", res[2][i])

        return bRunSucessful

    def exit_model(self):
        """ Exit the model (HEC-RAS specific)

        Returns
        -------

        """

        if self._project is not None:
            self._project = None

        if self._RASController is not None:
            print("Quitting HEC-RAS ...")
            self._RASController.QuitRas()

            del self._RASController

            print("Finished quitting HEC-RAS.")

    def get_plan_names(self, IncludeOnlyPlansInBaseDirectory):
        """Get the list of plan names in the current project

        Based on this COM object function:
        Plan_Names(self, PlanCount=defaultNamedNotOptArg, PlanNames=defaultNamedNotOptArg,
        IncludeOnlyPlansInBaseDirectory=defaultNamedNotOptArg)

        Parameters
        ----------
        IncludeOnlyPlansInBaseDirectory : bool
            whether only include plans in the base directory

        Returns
        -------
        PlanCount : int
            number of plans in the current project
        PlanNames : list
            a list of plan names in the current project

        """

        if self._RASController is None:
            print("RASController has not been created yet. Call init_model() first to create a RASController.")

        #it returns PlanCount, PlanNames and IncludeOnlyPlansInBaseDirectory (temp; not used)
        PlanCount, PlanNames, temp = self._RASController.Plan_Names(None, None, IncludeOnlyPlansInBaseDirectory)
        
        return PlanCount, PlanNames