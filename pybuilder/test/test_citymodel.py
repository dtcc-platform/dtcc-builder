import unittest
import sys
from pathlib import Path

sys.path.append( str((Path(__file__).parent / "..").resolve() ))

import _pybuilder
import Parameters
import PointCloud
import CityModelGenerator

data_dir = (Path(__file__).parent / "../../unittests/data").resolve()
p = Parameters.loadParameters()

class TestCityModel(unittest.TestCase):

    def setUp(self):
        shp_file = data_dir / "MinimalCase" / "PropertyMap.shp"
        self.cm = CityModelGenerator.generateCityModel(shp_file,p)
        self.pc = PointCloud.readLasFiles(data_dir / "MinimalCase" / "pointcloud.las")

    def test_bounds(self):
        shp_file = data_dir / "MinimalCase" / "PropertyMap.shp"
        bbox = CityModelGenerator.buildingBounds(shp_file)
        px, py, qx, qy = bbox
        self.assertAlmostEqual(px,-5.14247441879)
        self.assertAlmostEqual(qy, -1.09814696217)

    def test_padded_bounds(self):
        shp_file = data_dir / "MinimalCase" / "PropertyMap.shp"
        bbox = CityModelGenerator.buildingBounds(shp_file,10)
        px, py, qx, qy = bbox
        self.assertAlmostEqual(px,-15.14247441879)
        self.assertAlmostEqual(qy, 8.90185303782)
       

    def test_load(self):
        shp_file = data_dir / "MinimalCase" / "PropertyMap.shp"
        cm = CityModelGenerator.generateCityModel(shp_file,p)
        self.assertIsInstance(cm, _pybuilder.CityModel)
    
    def test_clean(self):
        shp_file = data_dir / "MinimalCase" / "PropertyMap.shp"
        cm = CityModelGenerator.generateCityModel(shp_file,p)
        clean_cm = CityModelGenerator.cleanCityModel(cm,0.5)

    def test_extarct_points(self):
        cm_with_points = CityModelGenerator.extractBuildingPoints(self.cm,self.pc,1.0,2.0)
        self.assertIsInstance(cm_with_points,_pybuilder.CityModel)

if __name__ == '__main__':
    unittest.main()