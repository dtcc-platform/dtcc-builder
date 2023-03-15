import unittest
import sys, os
from pathlib import Path

sys.path.append(str((Path(__file__).parent / "..").resolve()))

from dtcc_builder import _pybuilder

from dtcc_builder.ElevationModel import ElevationModel
from dtcc_builder.PointCloud import PointCloud
from dtcc_builder.CityModel import CityModel
from dtcc_builder.Parameters import load_parameters
from dtcc_builder.Utils import bounds_intersect
from dtcc_builder import Meshing


import dtcc_io as io


data_dir = (Path(__file__).parent / "../data").resolve()



class TestMeshing(unittest.TestCase):
    def setUp(self):
        shp_file = data_dir / "MinimalCase" / "PropertyMap.shp"
        self.p = load_parameters(data_dir / "MinimalCase" / "Parameters.json")


        self.cm_bounds = io.citymodel.building_bounds(shp_file, 10)
        self.pc_bounds = io.pointcloud.calc_las_bounds(data_dir / "MinimalCase")
        self.bounds = bounds_intersect(self.cm_bounds, self.pc_bounds)
        self.cm = CityModel(shp_file, parameters=self.p)
        self.cm.set_origin((self.bounds[0], self.bounds[1]))

        self.pc = PointCloud(data_dir / "MinimalCase" / "pointcloud.las")
        self.pc.set_origin((self.bounds[0], self.bounds[1]))

        self.cm.extract_building_points(self.pc, 1.0, 2.0)

        self.dtm = ElevationModel(self.pc, 0.5, [2, 9])
        self.cm.compute_building_heights(self.dtm, 0.9, 0.95)

        self.cm.clean_citymodel(self.p["MinVertexDistance"])

    def test_buildings_have_height(self):
        self.assertEqual(self.cm.buildings[0].height, 5)

    def test_mesh2D(self):
        mesh = Meshing.generate_mesh2D(self.cm, self.dtm.bounds, self.p["MeshResolution"])
        self.assertIsInstance(mesh, _pybuilder.Mesh2D)

    def test_mesh3D(self):
        mesh = Meshing.generate_mesh2D(self.cm, self.dtm.bounds, self.p["MeshResolution"])
        mesh3D = Meshing.generate_mesh3D(mesh, self.p["DomainHeight"], self.p["MeshResolution"])
        self.assertIsInstance(mesh3D, _pybuilder.Mesh3D)

    def test_smooth_mesh_no_fix(self):
        mesh = Meshing.generate_mesh2D(self.cm, self.dtm.bounds, self.p["MeshResolution"])
        mesh3D = Meshing.generate_mesh3D(mesh, self.p["DomainHeight"], self.p["MeshResolution"])
        smooth_mesh = Meshing.smooth_mesh3D(
            mesh3D, self.cm, self.dtm, self.p["DomainHeight"], False
        )
        self.assertIsInstance(smooth_mesh, _pybuilder.Mesh3D)

    def test_smooth_mesh_fix(self):
        mesh = Meshing.generate_mesh2D(self.cm, self.dtm.bounds, self.p["MeshResolution"])
        mesh3D = Meshing.generate_mesh3D(mesh, self.p["DomainHeight"], self.p["MeshResolution"])
        smooth_mesh = Meshing.smooth_mesh3D(
            mesh3D, self.cm, self.dtm, self.p["DomainHeight"], True
        )
        self.assertIsInstance(smooth_mesh, _pybuilder.Mesh3D)

    def test_generate_surface3D(self):
        surfaces = Meshing.generate_surface3D(self.cm, self.dtm, self.p["MeshResolution"])
        self.assertIsInstance(surfaces[0], _pybuilder.Surface3D)
        self.assertIsInstance(surfaces[-1], _pybuilder.Surface3D)

    def test_merge_surface3D(self):
        surfaces = Meshing.generate_surface3D(self.cm, self.dtm, self.p["MeshResolution"])
        merged_surface = Meshing.merge_surfaces3D(surfaces)
        self.assertIsInstance(merged_surface, _pybuilder.Surface3D)

    def test_trim_mesh(self):
        mesh2D = Meshing.generate_mesh2D(self.cm, self.dtm.bounds, self.p["MeshResolution"])
        mesh3D = Meshing.generate_mesh3D(mesh2D, self.p["DomainHeight"], self.p["MeshResolution"])

        trimmedMesh = Meshing.trim_mesh3D(mesh3D, mesh2D, self.cm, mesh3D.numLayers)
        self.assertIsInstance(trimmedMesh, _pybuilder.Mesh3D)

    def test_extract_boundary3D(self):
        mesh2D = Meshing.generate_mesh2D(self.cm, self.dtm.bounds, self.p["MeshResolution"])
        mesh3D = Meshing.generate_mesh3D(mesh2D, self.p["DomainHeight"], self.p["MeshResolution"])
        boundary = Meshing.extract_boundary3D(mesh3D)
        self.assertIsInstance(boundary, _pybuilder.Surface3D)

    def test_extract_open_surface(self):
        mesh2D = Meshing.generate_mesh2D(self.cm, self.dtm.bounds, self.p["MeshResolution"])
        mesh3D = Meshing.generate_mesh3D(mesh2D, self.p["DomainHeight"], self.p["MeshResolution"])
        boundary = Meshing.extract_boundary3D(mesh3D)
        open_surface = Meshing.extract_open_surface3D(boundary)
        self.assertIsInstance(open_surface, _pybuilder.Surface3D)

    # def test_write_mesh2D(self):
    #     mesh2D = Meshing.generate_mesh2D(self.cm, self.dtm.bounds, self.p["MeshResolution"])
    #     out_file = "test_2D.vtk"
    #     done = Meshing.write_VTK_mesh2D(mesh2D, out_file)
    #     self.assertTrue(done)
    #     self.assertTrue(os.path.isfile(out_file))
    #     os.unlink(out_file)

    def test_write_mesh3D(self):
        mesh2D = Meshing.generate_mesh2D(self.cm, self.dtm.bounds, self.p["MeshResolution"])
        mesh3D = Meshing.generate_mesh3D(mesh2D, self.p["DomainHeight"], self.p["MeshResolution"])

        out_file = "test_3D.vtk"
        done = Meshing.write_volume_mesh(mesh3D, out_file)
        self.assertTrue(os.path.isfile(out_file))
        os.unlink(out_file)

    def test_surface_mesh(self):
        mesh2D = Meshing.generate_mesh2D(self.cm, self.dtm.bounds, self.p["MeshResolution"])
        mesh3D = Meshing.generate_mesh3D(mesh2D, self.p["DomainHeight"], self.p["MeshResolution"])
        boundary = Meshing.extract_boundary3D(mesh3D)

        out_file = "test.obj"
        done = Meshing.write_surface(boundary, out_file)
        self.assertTrue(os.path.isfile(out_file))
        os.unlink(out_file)
