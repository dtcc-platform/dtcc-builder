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
        self.assertEqual(len(builder_pc.points), 8148)


class TestBuilderCityModel(unittest.TestCase):
    def test_convert_citymodel(self):
        cm = io.load_citymodel(data_dir / "MinimalCase" / "PropertyMap.shp")
        builder_cm = builder.builder_datamodel.create_builder_citymodel(cm)
        self.assertEqual(len(builder_cm.buildings), len(cm.buildings))
        self.assertEqual(len(builder_cm.buildings), 5)


class TestBuilderMesh(unittest.TestCase):
    def test_convert_mesh(self):
        mesh = io.load_mesh(data_dir / ""cube.vtk")
        builder_mesh = builder.builder_datamodel.create_builder_mesh(mesh)
        
        self.assertEqual(len(builder_mesh.vertices), len(mesh.vertices))
        self.assertEqual(len(builder_mesh.vertices), 8)
        self.assertEqual(len(builder_mesh.faces), len(mesh.faces))
        self.assertEqual(len(builder_mesh.faces), 12)
