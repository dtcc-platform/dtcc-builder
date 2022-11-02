import sys
from pathlib import Path
sys.path.append( str((Path(__file__).parent / "..").resolve() ))

import unittest

import _pybuilder
import Parameters
import PointCloud

data_dir = (Path(__file__).parent / "../../unittests/data").resolve()
p = Parameters.loadParameters()

class TestLoadPointCloud(unittest.TestCase):

    def test_load_dir(self):
        pc = PointCloud.readLasFiles(data_dir / "MinimalCase")
        self.assertIsInstance(pc, _pybuilder.PointCloud)
    def test_load_file(self):
        pc = PointCloud.readLasFiles(data_dir / "MinimalCase" / "pointcloud.las")
        self.assertIsInstance(pc, _pybuilder.PointCloud)
    def test_las_bounds(self):
        bbox = PointCloud.getLasBounds(data_dir / "MinimalCase")
        self.assertIsNotNone(bbox)
        px,py,qx,qy = bbox
        self.assertAlmostEqual(px,-8.017474418)
        self.assertAlmostEqual(qy,1.838260469410343)

class TestOutlierRemover(unittest.TestCase):
    def setUp(self) -> None:
        self.pc = PointCloud.readLasFiles(data_dir / "MinimalCase" / "pointcloud.las")
    
    def test_remove_outliers(self):
        filtered_pc = PointCloud.globalOutlierRemover(self.pc,1)
        self.assertIsInstance(filtered_pc,_pybuilder.PointCloud)
        
        
if __name__ == '__main__':
    unittest.main()