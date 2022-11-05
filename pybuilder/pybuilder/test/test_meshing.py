import unittest
import sys, os
from pathlib import Path

sys.path.append(str((Path(__file__).parent / "..").resolve()))

import _pybuilder
import ElevationModel
import PointCloud
import CityModel
import Parameters

import Meshing


data_dir = (Path(__file__).parent / "../../../unittests/data").resolve()
p = Parameters.loadParameters()

class TestMeshing(unittest.TestCase):
    def setUp(self):
        shp_file = data_dir / "MinimalCase" / "PropertyMap.shp"
        self.cm = CityModel.generateCityModel(shp_file, p)
        self.bounds = CityModel.buildingBounds(shp_file,p["DomainMargin"])
        self.pc = PointCloud.readLasFiles(data_dir / "MinimalCase" / "pointcloud.las")
        self.dtm = ElevationModel.generateElevationModel(self.pc, 0.5, [2, 9])

    def test_mesh2D(self):
        mesh = Meshing.generateMesh2D(self.cm, self.bounds, p["MeshResolution"] )
        self.assertIsInstance(mesh, _pybuilder.Mesh2D)