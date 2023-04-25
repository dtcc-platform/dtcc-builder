import unittest

import sys, os
from pathlib import Path
import dtcc_builder as builder
import dtcc_io as io

data_dir = (Path(__file__).parent / "../data").resolve()


class TestBuilderPointCloud(unittest.TestCase):
    def test_convert_pointcloud(self):
        pc = io.load_pointcloud(data_dir / "MinimalCase" / "pointcloud.las")
        builder_pc = builder.builder_datamodel.create_builder_pointcloud(pc)
        self.assertEqual(len(builder_pc.points), len(pc.points))
