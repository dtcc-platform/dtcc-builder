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
p = Parameters.loadParameters()

class TestMeshing(unittest.TestCase):
    def setUp(self):
        shp_file = data_dir / "MinimalCase" / "PropertyMap.shp"

        self.cm_bounds = CityModel.buildingBounds(shp_file,10)
        self.pc_bounds = PointCloud.getLasBounds(data_dir / "MinimalCase" )
        print(self.cm_bounds)
        print(self.pc_bounds)
        self.bounds = boundsIntersect(self.cm_bounds,self.pc_bounds)

        self.cm = CityModel.generateCityModel(shp_file, p)
        self.cm = CityModel.setOrigin(self.cm,(self.bounds[0],self.bounds[1]))
        
        self.pc = PointCloud.readLasFiles(data_dir / "MinimalCase" / "pointcloud.las")
        self.pc = PointCloud.setOrigin(self.pc,(self.bounds[0],self.bounds[1]))

        self.dtm = ElevationModel.generateElevationModel(self.pc, 0.5, [2, 9])



    def test_mesh2D(self):
        mesh = Meshing.generateMesh2D(self.cm, ElevationModel.getBounds(self.dtm), p["MeshResolution"] )
        self.assertIsInstance(mesh, _pybuilder.Mesh2D)