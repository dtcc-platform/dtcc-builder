import sys
from pathlib import Path

sys.path.append(str((Path(__file__).parent / "..").resolve()))

import unittest

import _pybuilder
import ElevationModel
import Parameters
import PointCloud

data_dir = (Path(__file__).parent / "../../../unittests/data").resolve()
p = Parameters.load_parameters()


class TestCreateDEM(unittest.TestCase):
    def setUp(self):
        self.pc = PointCloud.read_las_files(data_dir / "MinimalCase" / "pointcloud.las")
        self.dem = ElevationModel.generate_elevation_model(self.pc, 0.5, [2, 9])

    def test_generate_dem(self):
        self.assertIsInstance(self.dem, _pybuilder.GridField2D)

    def test_smooth_elevation(self):
        smooth_dem = ElevationModel.smooth_elevation_model(self.dem, 5)
        self.assertIsInstance(smooth_dem, _pybuilder.GridField2D)

    def test_get_bounds(self):
        bounds = ElevationModel.get_bounds(self.dem)
        self.assertEquals(
            bounds,
            (
                -8.01747441879072,
                -18.8503320990088,
                15.923729412288978,
                1.8382604694103435,
            ),
        )
