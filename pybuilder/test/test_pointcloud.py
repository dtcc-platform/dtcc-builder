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
        self.assertAlmostEqual(qy, 8.90185303782)

if __name__ == '__main__':
    unittest.main()