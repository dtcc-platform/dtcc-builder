import sys
from pathlib import Path

sys.path.append(str((Path(__file__).parent / "..").resolve()))

import unittest

from dtcc.builder import _pybuilder

from dtcc.builder import ElevationModel
from dtccpybuilder.Parameters import load_parameters
from dtccpybuilder.PointCloud import PointCloud

data_dir = (Path(__file__).parent / "../../../unittests/data").resolve()
p = load_parameters()


class TestCreateDEM(unittest.TestCase):
    def setUp(self):
        self.pc = PointCloud(data_dir / "MinimalCase" / "pointcloud.las")
        self.dem = ElevationModel(self.pc, 0.5, [2, 9])

    def test_generate_dem(self):
        self.assertIsInstance(self.dem._grid_field, _pybuilder.GridField2D)

    def test_smooth_elevation(self):
        smooth_dem = self.dem.smooth_elevation_model(5)
        self.assertIsInstance(self.dem._grid_field, _pybuilder.GridField2D)

    def test_get_bounds(self):
        bounds = self.dem.bounds

        self.assertAlmostEquals(bounds[0], -8.0174744, places=4)
        self.assertAlmostEquals(bounds[1], -18.850332, places=4)
        self.assertAlmostEquals(bounds[2], 15.923729412, places=4)
        self.assertAlmostEquals(bounds[3], 1.838260469, places=4)
