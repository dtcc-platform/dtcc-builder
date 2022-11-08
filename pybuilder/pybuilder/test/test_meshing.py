import unittest
import sys, os
from pathlib import Path

sys.path.append(str((Path(__file__).parent / "..").resolve()))

import _pybuilder
import ElevationModel
import PointCloud
import CityModel
import Parameters
from Utils import boundsIntersect
import Meshing


data_dir = (Path(__file__).parent / "../../../unittests/data").resolve()
p = Parameters.loadParameters(data_dir / "MinimalCase" / "Parameters.json")


class TestMeshing(unittest.TestCase):
    def setUp(self):
        shp_file = data_dir / "MinimalCase" / "PropertyMap.shp"

        self.cm_bounds = CityModel.buildingBounds(shp_file, 10)
        self.pc_bounds = PointCloud.getLasBounds(data_dir / "MinimalCase")
        print(self.cm_bounds)
        print(self.pc_bounds)
        self.bounds = boundsIntersect(self.cm_bounds, self.pc_bounds)

        self.cm = CityModel.generateCityModel(shp_file, p)
        self.cm = CityModel.setOrigin(self.cm, (self.bounds[0], self.bounds[1]))

        self.pc = PointCloud.readLasFiles(data_dir / "MinimalCase" / "pointcloud.las")
        self.pc = PointCloud.setOrigin(self.pc, (self.bounds[0], self.bounds[1]))

        self.cm = CityModel.extractBuildingPoints(self.cm, self.pc, 1.0, 2.0)

        self.dtm = ElevationModel.generateElevationModel(self.pc, 0.5, [2, 9])
        self.cm = CityModel.computeBuildingHeights(self.cm, self.dtm, 0.9, 0.95)

        self.cm = CityModel.cleanCityModel(self.cm, p["MinVertexDistance"])

    def test_buildings_have_height(self):
        self.assertEqual(self.cm.buildings[0].height, 5)

    def test_mesh2D(self):
        mesh = Meshing.generateMesh2D(
            self.cm, ElevationModel.getBounds(self.dtm), p["MeshResolution"]
        )
        self.assertIsInstance(mesh, _pybuilder.Mesh2D)

    def test_mesh3D(self):
        mesh = Meshing.generateMesh2D(
            self.cm, ElevationModel.getBounds(self.dtm), p["MeshResolution"]
        )
        mesh3D = Meshing.generateMesh3D(mesh, p["DomainHeight"], p["MeshResolution"])
        self.assertIsInstance(mesh3D, _pybuilder.Mesh3D)

    def test_smooth_mesh_no_fix(self):
        mesh = Meshing.generateMesh2D(
            self.cm, ElevationModel.getBounds(self.dtm), p["MeshResolution"]
        )
        mesh3D = Meshing.generateMesh3D(mesh, p["DomainHeight"], p["MeshResolution"])
        smooth_mesh = Meshing.smoothMesh3D(
            mesh3D, self.cm, self.dtm, p["DomainHeight"], False
        )
        self.assertIsInstance(smooth_mesh, _pybuilder.Mesh3D)

    def test_smooth_mesh_fix(self):
        mesh = Meshing.generateMesh2D(
            self.cm, ElevationModel.getBounds(self.dtm), p["MeshResolution"]
        )
        mesh3D = Meshing.generateMesh3D(mesh, p["DomainHeight"], p["MeshResolution"])
        smooth_mesh = Meshing.smoothMesh3D(
            mesh3D, self.cm, self.dtm, p["DomainHeight"], True
        )
        self.assertIsInstance(smooth_mesh, _pybuilder.Mesh3D)

    def test_generate_surface3D(self):
        surfaces = Meshing.generateSurface3D(self.cm, self.dtm, p["MeshResolution"])
        self.assertIsInstance(surfaces[0], _pybuilder.Surface3D)
        self.assertIsInstance(surfaces[-1], _pybuilder.Surface3D)

    def test_merge_surface3D(self):
        surfaces = Meshing.generateSurface3D(self.cm, self.dtm, p["MeshResolution"])
        merged_surface = Meshing.mergeSurfaces3D(surfaces)
        self.assertIsInstance(merged_surface, _pybuilder.Surface3D)

    def test_trim_mesh(self):
        mesh2D = Meshing.generateMesh2D(
            self.cm, ElevationModel.getBounds(self.dtm), p["MeshResolution"]
        )
        mesh3D = Meshing.generateMesh3D(mesh2D, p["DomainHeight"], p["MeshResolution"])

        trimmedMesh = Meshing.trimMesh3D(mesh3D, mesh2D, self.cm, mesh3D.numLayers)
        self.assertIsInstance(trimmedMesh, _pybuilder.Mesh3D)

    def test_extract_boundary3D(self):
        mesh2D = Meshing.generateMesh2D(
            self.cm, ElevationModel.getBounds(self.dtm), p["MeshResolution"]
        )
        mesh3D = Meshing.generateMesh3D(mesh2D, p["DomainHeight"], p["MeshResolution"])
        boundary = Meshing.extractBoundary3D(mesh3D)
        self.assertIsInstance(boundary, _pybuilder.Surface3D)

    def test_extract_open_surface(self):
        mesh2D = Meshing.generateMesh2D(
            self.cm, ElevationModel.getBounds(self.dtm), p["MeshResolution"]
        )
        mesh3D = Meshing.generateMesh3D(mesh2D, p["DomainHeight"], p["MeshResolution"])
        boundary = Meshing.extractBoundary3D(mesh3D)
        open_surface = Meshing.extractOpenSurface3D(boundary)
        self.assertIsInstance(open_surface, _pybuilder.Surface3D)
