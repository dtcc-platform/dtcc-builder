from dtcc_builder.meshing import (
    extrude_building,
    mesh_multisurface,
    mesh_surface,
    mesh_multisurfaces,
)
import dtcc_io as io
from dtcc_model import Mesh, Building, Surface, MultiSurface
from shapely.geometry import Polygon
import unittest
from pathlib import Path


project_dir = (Path(__file__).parent / "../data" / "MinimalCase").resolve()


class TestExtrudeBuilding(unittest.TestCase):
    def test_extrude_building(self):
        city = io.load_city(project_dir / "PropertyMap.shp")
        building = city.buildings[0]
        building.ground_level = 1
        building.height = 5
        mesh = extrude_building(building, max_mesh_size=10, min_mesh_angle=25)
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
        mesh = extrude_building(
            building, max_mesh_size=10, min_mesh_angle=25, ground_to_zero=True
        )
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
        mesh = extrude_building(
            building, max_mesh_size=10, min_mesh_angle=25, cap_base=True
        )
        self.assertIsInstance(mesh, Mesh)
        self.assertEqual(len(mesh.vertices), 12)
        self.assertEqual(len(mesh.faces), 12)

    def test_mesh_per_floor(self):
        building = Building(
            footprint=Polygon([(0, 0), (0, 10), (10, 10), (10, 0)]),
            height=20,
            floors=4,
            ground_level=5,
        )
        mesh = extrude_building(
            building, max_mesh_size=10, min_mesh_angle=25, per_floor=True, cap_base=True
        )
        self.assertEqual(
            len(mesh.faces), 8 * 4 + 5 * 2
        )  # 8 wall faces per floor, 2 roof faces per floor + 2 for the base
        self.assertEqual(
            len(mesh.vertices), 8 * 4 + 4
        )  # 8 vertices per floor + 4 for the base


class TestMeshSurfaces(unittest.TestCase):
    def test_mesh_simple_surface(self):
        surface = Surface(
            vertices=[
                [0, 0, 5],
                [0, 10, 5],
                [10, 10, 8],
                [10, 0, 8],
            ]
        )
        mesh = mesh_surface(surface)
        self.assertEqual(len(mesh.vertices), 4)
        self.assertEqual(len(mesh.faces), 2)
        self.assertEqual(mesh.vertices[:, 2].min(), 5)
        self.assertEqual(mesh.vertices[:, 2].max(), 8)

    def test_mesh_triangle_surface(self):
        surface = Surface(
            vertices=[
                [0, 0, 5],
                [0, 10, 5],
                [10, 10, 8],
                [10, 0, 8],
            ]
        )
        mesh = mesh_surface(surface, max_mesh_edge_size=5)
        self.assertEqual(len(mesh.vertices), 13)
        self.assertEqual(len(mesh.faces), 16)
        self.assertAlmostEqual(mesh.vertices[:, 2].min(), 5)
        self.assertAlmostEqual(mesh.vertices[:, 2].max(), 8)


class TestMeshMultiSurface(unittest.TestCase):
    def test_mesh_multisurface(self):
        surface1 = Surface(
            vertices=[
                [0, 0, 5],
                [0, 10, 5],
                [10, 10, 8],
                [10, 0, 8],
            ]
        )

        surface2 = Surface(
            vertices=[
                [0, 0, 0],
                [10, 0, 0],
                [10, 10, 0],
                [2, 10, 0],
                [2, 12, 0],
                [0, 12, 0],
            ]
        )
        multisurface = MultiSurface()
        multisurface.surfaces = [surface1, surface2]
        mesh = mesh_multisurface(multisurface)
        self.assertEqual(len(mesh.vertices), 10)
        self.assertEqual(len(mesh.faces), 6)
        self.assertAlmostEqual(mesh.vertices[:, 2].min(), 0)
        self.assertAlmostEqual(mesh.vertices[:, 2].max(), 8)


class TestMeshMultiSurfaces(unittest.TestCase):
    def test_mesh_multisurfaces(self):
        surface1 = Surface(
            vertices=[
                [0, 0, 5],
                [0, 10, 5],
                [10, 10, 8],
                [10, 0, 8],
            ]
        )

        surface2 = Surface(
            vertices=[
                [0, 0, 1],
                [10, 0, 1],
                [10, 10, 1],
                [2, 10, 1],
                [2, 12, 1],
                [0, 12, 1],
            ]
        )
        multisurface1 = MultiSurface()
        multisurface1.surfaces = [surface1, surface2]

        surface3 = Surface(
            vertices=[
                [0, 0, 4],
                [0, 10, 4],
                [10, 10, 7],
                [10, 0, 7],
            ]
        )

        surface4 = Surface(
            vertices=[
                [0, 0, 0],
                [10, 0, 0],
                [10, 10, 0],
                [2, 10, 0],
                [2, 12, 0],
                [0, 12, 0],
            ]
        )
        multisurface2 = MultiSurface()
        multisurface2.surfaces = [surface3, surface4]

        meshes = mesh_multisurfaces([multisurface1, multisurface2])
        self.assertEqual(len(meshes), 2)
        self.assertEqual(len(meshes[0].vertices), 10)
        self.assertEqual(len(meshes[0].faces), 6)
        self.assertAlmostEqual(meshes[0].vertices[:, 2].min(), 1)
        self.assertAlmostEqual(meshes[0].vertices[:, 2].max(), 8)
        self.assertEqual(len(meshes[1].vertices), 10)
        self.assertEqual(len(meshes[1].faces), 6)
        self.assertAlmostEqual(meshes[1].vertices[:, 2].min(), 0)
        self.assertAlmostEqual(meshes[1].vertices[:, 2].max(), 7)
