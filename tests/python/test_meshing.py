from dtcc_builder.meshing import extrude_building
import dtcc_io as io
from dtcc_model import Mesh

import unittest
from pathlib import Path


project_dir = (Path(__file__).parent / "../data" / "MinimalCase").resolve()


class TestExtrudeBuilding(unittest.TestCase):
    def test_extrude_building(self):
        city = io.load_city(project_dir / "PropertyMap.shp")
        building = city.buildings[0]
        building.ground_level = 1
        building.height = 5
        mesh = extrude_building(building, resolution=10)
        self.assertIsInstance(mesh, Mesh)
        self.assertEqual(len(mesh.vertices), 8)
        self.assertEqual(len(mesh.faces), 10)
        self.assertEqual(mesh.vertices[:, 2].min(), 1)
        self.assertEqual(mesh.vertices[:, 2].max(), 6)

    def test_extrude_building_ground_to_zero(self):
        city = io.load_city(project_dir / "PropertyMap.shp")
        building = city.buildings[0]
        building.ground_level = 1
        building.height = 5
        mesh = extrude_building(building, ground_to_zero=True)
        self.assertIsInstance(mesh, Mesh)
        self.assertEqual(len(mesh.vertices), 8)
        self.assertEqual(len(mesh.faces), 10)
        self.assertEqual(mesh.vertices[:, 2].min(), 0)
        self.assertEqual(mesh.vertices[:, 2].max(), 5)

    def test_extrude_building_capped(self):
        city = io.load_city(project_dir / "PropertyMap.shp")
        building = city.buildings[0]
        building.ground_level = 1
        building.height = 5
        mesh = extrude_building(building, cap_base=True)
        self.assertIsInstance(mesh, Mesh)
        self.assertEqual(len(mesh.vertices), 12)
        self.assertEqual(len(mesh.faces), 12)
