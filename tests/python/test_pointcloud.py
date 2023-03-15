import sys
from pathlib import Path

sys.path.append(str((Path(__file__).parent / "..").resolve()))

import unittest

import dtcc_io as io
from dtcc_builder import _pybuilder
from dtcc_builder import Parameters

from dtcc_builder import PointCloud

data_dir = (Path(__file__).parent / "../data").resolve()
p = Parameters.load_parameters(data_dir / "MinimalCase" / "Parameters.json")


class TestLoadPointCloud(unittest.TestCase):
    def test_load_dir(self):
        pc = PointCloud(data_dir / "MinimalCase")
        self.assertIsInstance(pc._builder_pc, _pybuilder.PointCloud)

    def test_load_file(self):
        pc = PointCloud(data_dir / "MinimalCase" / "pointcloud.las")
        self.assertIsInstance(pc._builder_pc, _pybuilder.PointCloud)
        self.assertEqual(len(pc), 8148)

    def test_load_file_bounded(self):
        pc = PointCloud(
            data_dir / "MinimalCase" / "pointcloud.las", bounds=(-2, -2, 0, 0)
        )
        self.assertIsInstance(pc._builder_pc, _pybuilder.PointCloud)
        self.assertEqual(len(pc), 64)

    def test_load_dir_bounded(self):
        pc = PointCloud(data_dir / "MinimalCase", bounds=(-2, -2, 0, 0))
        self.assertIsInstance(pc._builder_pc, _pybuilder.PointCloud)
        self.assertEqual(len(pc), 64)

    def test_las_bounds(self):
        bbox = io.pointcloud.calc_las_bounds(data_dir / "MinimalCase")
        self.assertIsNotNone(bbox)
        px, py, qx, qy = bbox
        self.assertAlmostEqual(px, -8.017474418)
        self.assertAlmostEqual(qy, 1.838260469)


class TestFilters(unittest.TestCase):
    def test_remove_outliers(self):
        pc = PointCloud(data_dir / "MinimalCase" / "pointcloud.las")
        pc.remove_global_outliers(1)
        self.assertIsInstance(pc._builder_pc, _pybuilder.PointCloud)
        self.assertEqual(len(pc), 8148 - 1183)

    def test_remove_vegetation(self):
        pc = PointCloud(data_dir / "MinimalCase" / "pointcloud.las")
        pc.vegetation_filter()
        self.assertIsInstance(pc._builder_pc, _pybuilder.PointCloud)


class TestConvert(unittest.TestCase):
    def test_convert2numpy(self):
        pc = PointCloud(data_dir / "MinimalCase" / "pointcloud.las")
        pts_array = pc.points
        l, d = pts_array.shape
        self.assertEqual(l, 8148)
        self.assertEqual(d, 3)

    def test_from_protobuf(self):
        pb_pc = str(data_dir / "MinimalCase" / "pointcloud.las.pb")
        with open(pb_pc, "rb") as src:
            pb_string = src.read()
        pc = PointCloud()
        pc.from_protobuf(pb_string)
        self.assertIsInstance(pc._builder_pc, _pybuilder.PointCloud)
        self.assertEquals(len(pc), 8148)

    def test_to_protobuf(self):
        pc = PointCloud(data_dir / "MinimalCase" / "pointcloud.las")
        protobuf = pc.to_protobuf()
        self.assertIsInstance(protobuf, bytes)

    def test_roundtrip_protobuf(self):
        pc = PointCloud(data_dir / "MinimalCase" / "pointcloud.las")
        protobuf = pc.to_protobuf()
        pc2 = PointCloud()
        pc2.from_protobuf(protobuf)
        self.assertIsInstance(pc2._builder_pc, _pybuilder.PointCloud)
        self.assertEquals(len(pc2), 8148)



if __name__ == "__main__":
    unittest.main()
