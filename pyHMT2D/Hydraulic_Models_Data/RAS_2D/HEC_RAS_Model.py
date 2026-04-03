"""
RAS_2D_Model:
"""

#import win32com.client as win32  #moved inside code
import os
import os.path

from pyHMT2D.Hydraulic_Models_Data import HydraulicModel

from pyHMT2D.__common__ import gVerbose

from pyHMT2D.Hydraulic_Models_Data.RAS_2D import RAS_2D_Data

from .helpers import *
import sys

import h5py


class HEC_RAS_RASMapper:
    """Parses a HEC-RAS .rasmap XML file to discover available terrain and land cover layers."""

    def __init__(self, rasmap_file):
        self.rasmap_file = rasmap_file
        self._project_dir = os.path.dirname(os.path.abspath(rasmap_file))
        self.terrains = []    # list of {"name": str, "filename": str (absolute path)}
        self.land_covers = [] # list of {"name": str, "filename": str (absolute path)}

    def parse(self):
        import xml.etree.ElementTree as ET
        tree = ET.parse(self.rasmap_file)
        root = tree.getroot()

        terrains_elem = root.find('Terrains')
        if terrains_elem is not None:
            for layer in terrains_elem.findall('Layer'):
                if layer.get('Type') == 'TerrainLayer':
                    filename = layer.get('Filename', '')
                    self.terrains.append({
                        'name': layer.get('Name', ''),
                        'filename': self._resolve_path(filename),
                    })

        map_layers_elem = root.find('MapLayers')
        if map_layers_elem is not None:
            for layer in map_layers_elem.findall('Layer'):
                if layer.get('Type') == 'LandCoverLayer':
                    filename = layer.get('Filename', '')
                    self.land_covers.append({
                        'name': layer.get('Name', ''),
                        'filename': self._resolve_path(filename),
                    })
        return self

    def get_first_terrain(self):
        """Return first terrain entry (document order), or None if no terrains defined."""
        return self.terrains[0] if self.terrains else None

    def _resolve_path(self, rel_path):
        """Resolve a rasmap relative path (.\\...) to an absolute path."""
        rel_path = rel_path.replace('\\', os.sep).replace('/', os.sep)
        if rel_path.startswith('.' + os.sep):
            rel_path = rel_path[2:]
        return os.path.normpath(os.path.join(self._project_dir, rel_path))

    def __repr__(self):
        return (f"HEC_RAS_RASMapper(terrains={[t['name'] for t in self.terrains]}, "
                f"land_covers={[lc['name'] for lc in self.land_covers]})")


class HEC_RAS_Geometry:
    """Data for a HEC-RAS geometry (.g## file and .g##.hdf file).

    The .g##.hdf file (geometry HDF) is created by HEC-RAS after geometry
    preprocessing and may not exist if preprocessing has not been run.
    This class reads the geometry HDF lazily (only when needed).
    """

    def __init__(self, geom_id, geom_file, project):
        """
        Parameters
        ----------
        geom_id : str
            Geometry short ID, e.g. 'g02'.
        geom_file : str
            Absolute path to the .g## text file.
        project : HEC_RAS_Project
            Back-reference to the owning project.
        """
        self.geom_id = geom_id
        self.geom_file = geom_file
        self.hdf_file = geom_file + '.hdf'
        self.project = project
        self._title = None
        self._terrain_hdf_file = None
        self._terrain_tiff_file = None
        self._terrain_resolved = False

    @property
    def title(self):
        if self._title is None:
            self._parse_geom_file()
        return self._title

    def _parse_geom_file(self):
        self._title = ''
        if not os.path.isfile(self.geom_file):
            return
        with open(self.geom_file, 'r', errors='replace') as f:
            for line in f:
                if line.startswith('Geom Title='):
                    self._title = line.split('=', 1)[1].strip()
                    break

    def hdf_exists(self):
        """Return True if the geometry HDF file exists."""
        return os.path.isfile(self.hdf_file)

    @property
    def terrain_tiff_file(self):
        """Absolute path to the terrain GeoTIFF file, resolved via fallback chain."""
        if not self._terrain_resolved:
            self._resolve_terrain()
        return self._terrain_tiff_file

    @property
    def terrain_hdf_file(self):
        """Absolute path to the RASMapper terrain HDF file."""
        if not self._terrain_resolved:
            self._resolve_terrain()
        return self._terrain_hdf_file

    def _resolve_terrain(self):
        """Resolve terrain files using the fallback chain:
        1. Geometry HDF Terrain Filename attribute
        2. First terrain in project .rasmap (document order)
        3. First .tif file found in the Terrain/ subdirectory
        4. None (warn; face-point Z will be NaN)
        """
        self._terrain_resolved = True

        # Step 1: geometry HDF
        if self.hdf_exists():
            try:
                hdf_path, tiff_path = extract_terrain_file_names(self.geom_file)
                self._terrain_hdf_file = hdf_path
                if tiff_path and os.path.isfile(tiff_path):
                    self._terrain_tiff_file = tiff_path
                    return
                elif tiff_path:
                    print(f"Warning: Terrain TIFF referenced in geometry HDF not found: {tiff_path}")
            except Exception as e:
                print(f"Warning: Could not read terrain from geometry HDF ({self.hdf_file}): {e}")

        # Step 2: first terrain in rasmap
        if self.project and self.project.rasmap:
            first = self.project.rasmap.get_first_terrain()
            if first:
                terrain_hdf = first['filename']
                self._terrain_hdf_file = terrain_hdf
                if os.path.isfile(terrain_hdf):
                    try:
                        tiff_name = get_terrain_file_attribute(terrain_hdf)
                        if tiff_name:
                            terrain_dir = os.path.dirname(terrain_hdf)
                            tiff_path = os.path.normpath(os.path.join(terrain_dir, tiff_name))
                            if os.path.isfile(tiff_path):
                                self._terrain_tiff_file = tiff_path
                                print(f"Info: Using terrain from rasmap for geometry {self.geom_id}: {tiff_path}")
                                return
                    except Exception as e:
                        print(f"Warning: Could not read terrain from rasmap entry ({terrain_hdf}): {e}")

        # Step 3: scan Terrain/ subdirectory
        project_dir = os.path.dirname(os.path.abspath(self.geom_file))
        terrain_dir = os.path.join(project_dir, 'Terrain')
        if os.path.isdir(terrain_dir):
            for fname in os.listdir(terrain_dir):
                if fname.lower().endswith('.tif') or fname.lower().endswith('.tiff'):
                    self._terrain_tiff_file = os.path.join(terrain_dir, fname)
                    print(f"Warning: Auto-discovered terrain TIFF for geometry {self.geom_id}: "
                          f"{self._terrain_tiff_file}. Verify this is the correct terrain file.")
                    return

        # Step 4: not found
        print(f"Warning: No terrain TIFF file found for geometry {self.geom_id}. "
              f"Face-point elevations will be set to NaN.")
        self._terrain_tiff_file = None

    def __str__(self):
        return f"{self.geom_id}: {self.title}"

    def __repr__(self):
        return (f"HEC_RAS_Geometry(geom_id='{self.geom_id}', "
                f"title='{self.title}', "
                f"hdf_exists={self.hdf_exists()})")


class HEC_RAS_SteadyFlow:
    """Data for a HEC-RAS steady flow (.f## file)."""

    def __init__(self, flow_id, flow_file, project):
        self.flow_id = flow_id
        self.flow_file = flow_file
        self.project = project
        self._title = None

    @property
    def title(self):
        if self._title is None:
            self._parse_flow_file()
        return self._title

    def _parse_flow_file(self):
        self._title = ''
        if not os.path.isfile(self.flow_file):
            return
        with open(self.flow_file, 'r', errors='replace') as f:
            for line in f:
                if line.startswith('Flow Title='):
                    self._title = line.split('=', 1)[1].strip()
                    break

    def __str__(self):
        return f"{self.flow_id}: {self.title}"

    def __repr__(self):
        return f"HEC_RAS_SteadyFlow(flow_id='{self.flow_id}', title='{self.title}')"


class HEC_RAS_UnsteadyFlow:
    """Data for a HEC-RAS unsteady flow (.u## file).

    Note: DSS (HEC Data Storage System) support for reading/writing boundary
    condition time series is a planned future feature.
    """

    def __init__(self, flow_id, flow_file, project):
        self.flow_id = flow_id
        self.flow_file = flow_file
        self.project = project
        self._title = None
        # TODO (future): DSS support for reading/writing boundary condition time series.

    @property
    def title(self):
        if self._title is None:
            self._parse_flow_file()
        return self._title

    def _parse_flow_file(self):
        self._title = ''
        if not os.path.isfile(self.flow_file):
            return
        with open(self.flow_file, 'r', errors='replace') as f:
            for line in f:
                if line.startswith('Flow Title='):
                    self._title = line.split('=', 1)[1].strip()
                    break

    def __str__(self):
        return f"{self.flow_id}: {self.title}"

    def __repr__(self):
        return f"HEC_RAS_UnsteadyFlow(flow_id='{self.flow_id}', title='{self.title}')"


class HEC_RAS_Plan:
    """Data for a HEC-RAS plan (.p## file).

    Each plan owns a RAS_2D_Data object (populated lazily via load_results()).
    """

    def __init__(self, plan_id, plan_file, project):
        """
        Parameters
        ----------
        plan_id : str
            Plan short ID, e.g. 'p03'.
        plan_file : str
            Absolute path to the .p## text file.
        project : HEC_RAS_Project
            Back-reference to the owning project.
        """
        self.plan_id = plan_id
        self.plan_file = plan_file
        self.hdf_file = plan_file + '.hdf'
        self.project = project
        self._plan_name = None
        self._steady = None
        self._geometry_id = None
        self._flow_id = None
        self.ras_2d_data = None  # populated by load_results()

    @property
    def plan_name(self):
        if self._plan_name is None:
            self._parse_plan_file()
        return self._plan_name

    @property
    def steady(self):
        if self._steady is None:
            self._parse_plan_file()
        return self._steady

    @property
    def geometry_id(self):
        if self._geometry_id is None:
            self._parse_plan_file()
        return self._geometry_id

    @property
    def flow_id(self):
        if self._flow_id is None:
            self._parse_plan_file()
        return self._flow_id

    @property
    def geometry(self):
        """Return the HEC_RAS_Geometry object associated with this plan."""
        if self.project and self.geometry_id:
            return self.project.geometries.get(self.geometry_id)
        return None

    def _parse_plan_file(self):
        self._plan_name = ''
        self._steady = False
        self._geometry_id = ''
        self._flow_id = ''
        if not os.path.isfile(self.plan_file):
            return
        has_unsteady = False
        with open(self.plan_file, 'r', errors='replace') as f:
            for line in f:
                line = line.strip()
                if line.startswith('Plan Title='):
                    self._plan_name = line.split('=', 1)[1].strip()
                elif line.startswith('Geom File='):
                    self._geometry_id = line.split('=', 1)[1].strip()
                elif line.startswith('Flow File='):
                    self._flow_id = line.split('=', 1)[1].strip()
                    fid = self._flow_id
                    if fid.lower().startswith('u'):
                        has_unsteady = True
        self._steady = not has_unsteady

    def hdf_exists(self):
        """Return True if the plan result HDF file exists."""
        return os.path.isfile(self.hdf_file)

    def load_results(self, terrain_file=None):
        """Load plan simulation results into ras_2d_data.

        Parameters
        ----------
        terrain_file : str, optional
            Path to a GeoTIFF terrain file. Overrides auto-detection from the
            geometry HDF and rasmap. Required only if auto-detection fails.

        Returns
        -------
        RAS_2D_Data
            The loaded data object (also stored as self.ras_2d_data).
        """
        from pyHMT2D.Hydraulic_Models_Data.RAS_2D.RAS_2D_Data import RAS_2D_Data
        if not self.hdf_exists():
            raise FileNotFoundError(
                f"Result HDF file '{self.hdf_file}' not found. "
                "Run the HEC-RAS simulation first."
            )
        self.ras_2d_data = RAS_2D_Data(self, terrain_file)
        return self.ras_2d_data

    def __str__(self):
        return f"{self.plan_id}: {self.plan_name}"

    def __repr__(self):
        return (f"HEC_RAS_Plan(plan_id='{self.plan_id}', "
                f"plan_name='{self.plan_name}', "
                f"geometry_id='{self.geometry_id}', "
                f"flow_id='{self.flow_id}', "
                f"steady={self.steady})")


class HEC_RAS_Project:
    """HEC-RAS project class — parses .prj and .rasmap files.

    Can be used standalone (no COM / HEC-RAS installation required) for
    post-processing simulation results.

    Typical usage (post-processing only)::

        project = HEC_RAS_Project("Muncie.prj")
        plan    = project.get_plan("p03")      # by ID or title
        data    = plan.load_results()           # terrain auto-detected
        data.saveHEC_RAS2D_results_to_VTK(...)

    When running simulations use HEC_RAS_Model (which extends this via
    composition and adds COM control).
    """

    def __init__(self, prj_file):
        """
        Parameters
        ----------
        prj_file : str
            Path to the HEC-RAS project file (.prj).
        """
        if not os.path.isfile(prj_file):
            raise FileNotFoundError(f"Project file '{prj_file}' not found.")

        self.prj_file = os.path.abspath(prj_file)
        self.project_dir = os.path.dirname(self.prj_file)
        self.project_name = os.path.splitext(os.path.basename(self.prj_file))[0]

        self.title = ''
        self.units = ''  # 'Feet' or 'Meter'

        # Ordered dicts: short-ID (str) -> object, e.g. "p03" -> HEC_RAS_Plan
        self.plans = {}
        self.geometries = {}
        self.unsteady_flows = {}
        self.steady_flows = {}

        self.current_plan_id = None  # e.g. "p04"

        # RASMapper layer pool (None if .rasmap file not found)
        self.rasmap = None

        self._parse_prj()
        self._parse_rasmap()

    def _parse_prj(self):
        geom_ids = []
        plan_ids = []
        unsteady_ids = []
        steady_ids = []

        with open(self.prj_file, 'r', errors='replace') as f:
            for line in f:
                line = line.strip()
                if line.startswith('Proj Title='):
                    self.title = line.split('=', 1)[1].strip()
                elif line == 'English Units':
                    self.units = 'Feet'
                elif line == 'SI Units':
                    self.units = 'Meter'
                elif line.startswith('Geom File='):
                    geom_ids.append(line.split('=', 1)[1].strip())
                elif line.startswith('Plan File='):
                    plan_ids.append(line.split('=', 1)[1].strip())
                elif line.startswith('Unsteady File='):
                    unsteady_ids.append(line.split('=', 1)[1].strip())
                elif line.startswith('Flow File='):
                    steady_ids.append(line.split('=', 1)[1].strip())
                elif line.startswith('Current Plan='):
                    self.current_plan_id = line.split('=', 1)[1].strip()

        for gid in geom_ids:
            geom_file = os.path.join(self.project_dir, self.project_name + '.' + gid)
            self.geometries[gid] = HEC_RAS_Geometry(gid, geom_file, self)

        for pid in plan_ids:
            plan_file = os.path.join(self.project_dir, self.project_name + '.' + pid)
            self.plans[pid] = HEC_RAS_Plan(pid, plan_file, self)

        for uid in unsteady_ids:
            flow_file = os.path.join(self.project_dir, self.project_name + '.' + uid)
            self.unsteady_flows[uid] = HEC_RAS_UnsteadyFlow(uid, flow_file, self)

        for fid in steady_ids:
            flow_file = os.path.join(self.project_dir, self.project_name + '.' + fid)
            self.steady_flows[fid] = HEC_RAS_SteadyFlow(fid, flow_file, self)

    def _parse_rasmap(self):
        rasmap_file = os.path.join(self.project_dir, self.project_name + '.rasmap')
        if os.path.isfile(rasmap_file):
            self.rasmap = HEC_RAS_RASMapper(rasmap_file)
            self.rasmap.parse()

    def get_plan(self, plan_id_or_name):
        """Get a plan by short ID (e.g. 'p03') or by plan title.

        Parameters
        ----------
        plan_id_or_name : str
            Plan short ID or plan title string.

        Returns
        -------
        HEC_RAS_Plan
        """
        if plan_id_or_name in self.plans:
            return self.plans[plan_id_or_name]
        for plan in self.plans.values():
            if plan.plan_name == plan_id_or_name:
                return plan
        raise KeyError(
            f"Plan '{plan_id_or_name}' not found. "
            f"Available plan IDs: {list(self.plans.keys())}"
        )

    @property
    def current_plan(self):
        """Return the current plan as an HEC_RAS_Plan object."""
        if self.current_plan_id and self.current_plan_id in self.plans:
            return self.plans[self.current_plan_id]
        return None

    def set_current_plan_id(self, plan_id):
        """Set current plan by short ID."""
        if plan_id not in self.plans:
            raise KeyError(f"Plan ID '{plan_id}' not found.")
        self.current_plan_id = plan_id

    def __str__(self):
        return self.title

    def __repr__(self):
        return (f"HEC_RAS_Project(title='{self.title}', "
                f"units='{self.units}', "
                f"plans={list(self.plans.keys())}, "
                f"geometries={list(self.geometries.keys())})")


class HEC_RAS_Model(HydraulicModel):
    """HEC-RAS Model class

    HEC-RAS Model controls the run of HEC-RAS. A plan can be loaded and executed.
    Currently HEC_RAS_Model has very limited capability to modify HEC-RAS project,
    plan, geometry, and flow data.     It is mainly used for Monte Carlo or batch simulations.

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
            if gVerbose:
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

        if gVerbose:
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
        elif (self.getVersion() == '6.6'):
            self._RASController = win32.gencache.EnsureDispatch('RAS66.HECRASController')
        else:
            raise Exception("The specified version of HEC-RAS is not currently supported.")

        # show HEC-RAS main window
        if not self._faceless:
            self._RASController.ShowRas()

        # print and confirm the HEC-RAS version
        if gVerbose:
            print('Successfully launched HEC-RAS version ', self._RASController.HECRASVersion())

        #_ = input("Press ENTER to continue:")


    def open_project(self, projectFileName):
        """Open the specified HEC-RAS project.

        Parses the project file structure and opens the project via the COM
        interface. The returned HEC_RAS_Project can also be created without
        COM for post-processing-only workflows.

        Parameters
        ----------
        projectFileName : str
            Path to the HEC-RAS project file (.prj).

        Returns
        -------
        HEC_RAS_Project
            The parsed project object.
        """
        if gVerbose:
            print("HEC-RAS opens project: ", projectFileName)

        if self._RASController is None:
            self.init_model()

        if os.path.isfile(projectFileName):
            self._project_file_name = os.path.abspath(projectFileName)
        else:
            raise FileNotFoundError(f'Project file "{projectFileName}" not found.')

        # Open via COM
        self._RASController.Project_Open(self._project_file_name)

        # Parse the project file structure (pure Python, no COM)
        self._project = HEC_RAS_Project(self._project_file_name)

        # Reconcile current plan with what COM reports
        current_plan_file = self._RASController.CurrentPlanFile()
        if current_plan_file:
            current_plan_file_abs = os.path.abspath(current_plan_file)
            for pid, plan in self._project.plans.items():
                if os.path.abspath(plan.plan_file) == current_plan_file_abs:
                    self._project.current_plan_id = pid
                    break

        if gVerbose:
            print("Project opened:", self._project)
            print("Current plan:", self._project.current_plan)

        return self._project

    def set_current_plan(self, plan_id_or_name):
        """Set the current plan by short ID (e.g. 'p03') or plan title.

        Parameters
        ----------
        plan_id_or_name : str
            Plan short ID (e.g. 'p03') or plan title string.
        """
        plan = self._project.get_plan(plan_id_or_name)
        # COM needs the plan title (name), not the ID
        ret = self._RASController.Plan_SetCurrent(plan.plan_name)
        if ret[0]:
            self._project.current_plan_id = plan.plan_id
            if gVerbose:
                print(f"Current plan set to: {plan}")
        else:
            print(f"Warning: COM Plan_SetCurrent failed for '{plan.plan_name}'. "
                  f"Updating internal state only.")
            self._project.current_plan_id = plan.plan_id

    def load_current_plan_results(self, terrain_file=None):
        """Load simulation results for the current plan.

        Parameters
        ----------
        terrain_file : str, optional
            Path to a GeoTIFF terrain file. Passed to plan.load_results().
            If None, terrain is auto-detected from the geometry HDF and rasmap.

        Returns
        -------
        RAS_2D_Data
        """
        plan = self._project.current_plan
        if plan is None:
            raise RuntimeError("No current plan set. Call open_project() first.")

        if gVerbose:
            print(f"Loading results for plan: {plan}")

        return plan.load_results(terrain_file=terrain_file)

    def get_simulation_case(self, bReload=False):
        """Return the RAS_2D_Data for the current plan, loading if necessary.

        Parameters
        ----------
        bReload : bool, optional
            If True, reload results from the HDF file even if already loaded.

        Returns
        -------
        RAS_2D_Data
        """
        plan = self._project.current_plan
        if plan is None:
            raise RuntimeError("No current plan set. Call open_project() first.")

        if plan.ras_2d_data is None or bReload:
            self.load_current_plan_results()

        return plan.ras_2d_data

    def get_current_planFile(self):
        """Return the current plan file path."""
        return self._RASController.CurrentPlanFile()

    def get_current_project(self):
        """Return the current HEC_RAS_Project."""
        if self._project is not None:
            return self._project
        raise RuntimeError("No project opened. Call open_project() first.")

    def save_project(self):
        """Save the current project

        Returns
        -------

        """

        if (self._RASController is not None)  and (self._project is not None):
            if gVerbose:
                print("Saving project: ", self._RASController.CurrentProjectTitle())
            self._RASController.Project_Save()

    def close_project(self):
        """Close the current project (if any)

        Returns
        -------

        """

        if (self._RASController is not None) and (self._project is not None):
            if gVerbose:
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

        if gVerbose:
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
            if gVerbose:
                print("HEC-RAS computed successfully.")
        else:
            if gVerbose:
                print("HEC-RAS computed unsuccessfully. The HEC-RAS Controller's Compute_CurrentPlan() function returned "
                      "False.")
            bRunSucessful = False

        if gVerbose:
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
            if gVerbose:
                print("Quitting HEC-RAS ...")
            self._RASController.QuitRas()

            del self._RASController

            if gVerbose:
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
