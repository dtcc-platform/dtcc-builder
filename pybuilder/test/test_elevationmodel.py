import sys
from pathlib import Path
sys.path.append( str((Path(__file__).parent / "..").resolve() ))

import unittest

import _pybuilder
import ElevationModel
import Parameters
import PointCloud

data_dir = (Path(__file__).parent / "../../unittests/data").resolve()
p = Parameters.loadParameters()

class TestCreateDEM(unittest.TestCase):
    def setUp(self):
        self.pc = PointCloud.readLasFiles(data_dir / "MinimalCase" / "pointcloud.las")
    
    def test_generate_dem(self):
        dem = ElevationModel.generateElevationModel(self.pc,0.5,[2,9])
        self.assertIsInstance(dem,_pybuilder.GridField2D)




