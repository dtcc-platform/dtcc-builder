import sys
from pathlib import Path

sys.path.append(str((Path(__file__).parent / "..").resolve()))

import unittest

import _pybuilder
import ElevationModel
import Parameters
import PointCloud

data_dir = (Path(__file__).parent / "../../../unittests/data").resolve()
p = Parameters.loadParameters()


class TestCreateDEM(unittest.TestCase):
    def setUp(self):
        self.pc = PointCloud.readLasFiles(data_dir / "MinimalCase" / "pointcloud.las")
        self.dem = ElevationModel.generateElevationModel(self.pc, 0.5, [2, 9])

    def test_generate_dem(self):
        self.assertIsInstance(self.dem, _pybuilder.GridField2D)

    def test_smooth_elevation(self):
        smooth_dem = ElevationModel.smoothElevationModel(self.dem,5)
        self.assertIsInstance(smooth_dem, _pybuilder.GridField2D)
    