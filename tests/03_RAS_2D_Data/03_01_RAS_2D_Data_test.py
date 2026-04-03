import pyHMT2D


class _MinimalGeom:
    def __init__(self, terrain_tiff_file):
        self.terrain_tiff_file = terrain_tiff_file
        self.terrain_hdf_file = None


class _MinimalPlan:
    """Minimal plan-like object for testing RAS_2D_Data without a full project."""
    def __init__(self, hdf_file, terrain_tiff_file=None):
        self.hdf_file = hdf_file
        self.geometry = None  # will trigger fallback
        self._terrain_tiff_file = terrain_tiff_file


def test_RAS_2D_Data_construction():
    """Test the construction of RAS_2D_Data object from a minimal plan object."""

    # pytest is run from the project root, so paths are relative to that.
    geom = _MinimalGeom("tests/data/TerrainMuncie_composite.tif")
    plan = _MinimalPlan("tests/data/Muncie2D.p01.hdf")
    plan.geometry = geom
    my_ras_2d_data = pyHMT2D.RAS_2D.RAS_2D_Data(plan)

    #print(my_ras_2d_data.TwoDAreaFace_FacePoints[0][0] )

    assert my_ras_2d_data.TwoDAreaFace_FacePoints[0][0][0] == 3
    assert my_ras_2d_data.TwoDAreaFace_FacePoints[0][0][1] == 2
