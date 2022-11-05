import sys
from pathlib import Path

sys.path.append(str((Path(__file__).parent / "..").resolve()))

import unittest

import _pybuilder
import Parameters
import PointCloud

data_dir = (Path(__file__).parent / "../../../unittests/data").resolve()
p = Parameters.loadParameters()


class TestLoadPointCloud(unittest.TestCase):
    def test_load_dir(self):
        pc = PointCloud.readLasFiles(data_dir / "MinimalCase")
        self.assertIsInstance(pc, _pybuilder.PointCloud)

    def test_load_file(self):
        pc = PointCloud.readLasFiles(data_dir / "MinimalCase" / "pointcloud.las")
        self.assertIsInstance(pc, _pybuilder.PointCloud)
        self.assertEqual(len(pc), 8148)

    def test_las_bounds(self):
        bbox = PointCloud.getLasBounds(data_dir / "MinimalCase")
        self.assertIsNotNone(bbox)
        px, py, qx, qy = bbox
        self.assertAlmostEqual(px, -8.017474418)
        self.assertAlmostEqual(qy, 1.838260469410343)


class TestFilters(unittest.TestCase):
    def test_remove_outliers(self):
        pc = PointCloud.readLasFiles(data_dir / "MinimalCase" / "pointcloud.las")
        # print(len(pc))
        filtered_pc = PointCloud.globalOutlierRemover(pc, 1)
        self.assertIsInstance(filtered_pc, _pybuilder.PointCloud)
        self.assertEqual(len(filtered_pc), 8148 - 1183)

    def test_remove_vegetation(self):
        pc = PointCloud.readLasFiles(data_dir / "MinimalCase" / "pointcloud.las")
        filtered_pc = PointCloud.VegetationFilter(pc)
        self.assertIsInstance(filtered_pc, _pybuilder.PointCloud)


class TestConvert(unittest.TestCase):
    def test_convert2numpy(self):
        pc = PointCloud.readLasFiles(data_dir / "MinimalCase" / "pointcloud.las")
        pts_array = PointCloud.pointCloud2numpy(pc)
        l, d = pts_array.shape
        self.assertEqual(l, 8148)
        self.assertEqual(d, 3)


if __name__ == "__main__":
    unittest.main()
