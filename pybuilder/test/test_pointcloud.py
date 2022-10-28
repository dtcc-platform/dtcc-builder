import sys
from pathlib import Path
sys.path.append( str((Path(__file__).parent / "..").resolve() ))

import unittest

import _pybuilder
import Parameters
import PointCloud

data_dir = (Path(__file__).parent / "../../unittests/data").resolve()
p = Parameters.loadParameters()

class TestPointCloud(unittest.TestCase):

    def test_load_dir(self):
        pc = PointCloud.ReadLasFiles(data_dir / "MinimalCase")
        self.assertIsInstance(pc, _pybuilder.PointCloud)
    def test_load_file(self):
        pc = PointCloud.ReadLasFiles(data_dir / "MinimalCase" / "pointcloud.las")
        self.assertIsInstance(pc, _pybuilder.PointCloud)
    # def test_las_bounds(self):
    #     bbox = PointCloud.GetLasBounds(data_dir / "MinimalCase")
    

if __name__ == '__main__':
    unittest.main()