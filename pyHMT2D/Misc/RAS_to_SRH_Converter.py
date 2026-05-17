"""
A Python class to convert HEC-RAS case, mainly mesh and material (Manning's n), to SRH-2D case
"""

import pyHMT2D

class RAS_to_SRH_Converter:
    """Converts a HEC-RAS 2D plan result to SRH-2D input format.

    For cases with multiple 2D areas, specify which area to convert via
    twoDAreaNumber (default: 0, the first area).
    """

    def __init__(self, prj_file, plan_id_or_name, srh_case_name,
                 twoDAreaNumber=0, terrain_file=None):
        """
        Parameters
        ----------
        prj_file : str
            Path to the HEC-RAS project file (.prj).
        plan_id_or_name : str
            Plan short ID (e.g. 'p03') or plan title to convert.
        srh_case_name : str
            Base name for the output SRH-2D files (.srhgeom, .srhmat).
        twoDAreaNumber : int, optional
            Index of the 2D area to convert (default 0, i.e. first area).
            For single-area cases this should always be 0.
        terrain_file : str, optional
            Path to a single-source GeoTIFF terrain file. Overrides
            auto-discovery from the geometry HDF / rasmap. Required when
            the project uses a composite terrain (stack of TIFFs), since
            auto-discovery currently picks only the first layer in the
            stack and produces incorrect floodplain elevations.
            See: https://github.com/psu-efd/pyHMT2D/issues/13
        """
        self.prj_file = prj_file
        self.plan_id_or_name = plan_id_or_name
        self.srh_case_name = srh_case_name
        self.twoDAreaNumber = twoDAreaNumber
        self.terrain_file = terrain_file

    def convert_to_SRH(self):
        """Perform the conversion and write .srhgeom and .srhmat files."""
        import pyHMT2D
        from pyHMT2D.Hydraulic_Models_Data.RAS_2D.HEC_RAS_Model import HEC_RAS_Project

        project = HEC_RAS_Project(self.prj_file)
        plan = project.get_plan(self.plan_id_or_name)

        if not plan.hdf_exists():
            raise FileNotFoundError(
                f"Plan result HDF '{plan.hdf_file}' not found. "
                "Run the HEC-RAS simulation first."
            )

        ras_2d_data = plan.load_results(terrain_file=self.terrain_file)

        n_areas = len(ras_2d_data.TwoDAreaNames)
        if self.twoDAreaNumber >= n_areas:
            raise IndexError(
                f"twoDAreaNumber={self.twoDAreaNumber} is out of range. "
                f"This plan has {n_areas} 2D area(s) (indices 0..{n_areas-1})."
            )

        # exportSRH{GEOM,MAT}File append the extension themselves — pass the
        # base name only to avoid producing files like "Muncie.srhgeom.srhgeom".
        ras_2d_data.exportSRHGEOMFile(self.srh_case_name,
                                       twoDAreaNumber=self.twoDAreaNumber)
        ras_2d_data.exportSRHMATFile(self.srh_case_name,
                                     twoDAreaNumber=self.twoDAreaNumber)
