import unittest
import sys, os
from pathlib import Path

sys.path.append(str((Path(__file__).parent / "..").resolve()))

from pybuilder import _pybuilder

from pybuilder.Parameters import load_parameters
from pybuilder.PointCloud import PointCloud
from pybuilder.ElevationModel import ElevationModel
from pybuilder.CityModel import CityModel, building_bounds, get_building_roof_points

data_dir = (Path(__file__).parent / "../../../unittests/data").resolve()
p = load_parameters(data_dir / "MinimalCase" / "Parameters.json")


class TestCityModel(unittest.TestCase):
    def setUp(self):
        self.shp_file = data_dir / "MinimalCase" / "PropertyMap.shp"
        self.cm = CityModel(self.shp_file, p)

        self.pc = PointCloud(data_dir / "MinimalCase" / "pointcloud.las")

    def test_bounds(self):
        bbox = building_bounds(self.shp_file)
        px, py, qx, qy = bbox
        self.assertAlmostEqual(px, -5.14247441879)
        self.assertAlmostEqual(qy, -1.09814696217)

    def test_padded_bounds(self):
        bbox = building_bounds(self.shp_file, 10)
        px, py, qx, qy = bbox
        self.assertAlmostEqual(px, -15.14247441879)
        self.assertAlmostEqual(qy, 8.90185303782)

    def test_set_origin(self):
        bbox = building_bounds(self.shp_file, 10)
        origin = (bbox[0], bbox[1])
        self.cm.set_origin(origin)
        self.assertEquals(self.cm.origin, origin)

    def test_load(self):
        shp_file = data_dir / "MinimalCase" / "PropertyMap.shp"
        cm = CityModel(shp_file, p)
        self.assertIsInstance(cm._builder_cm, _pybuilder.CityModel)

    def test_clean(self):

        self.cm.clean_citymodel(0.5)
        self.assertIsInstance(self.cm._builder_cm, _pybuilder.CityModel)
        self.assertTrue(self.cm.cleaned)

    def test_extarct_points(self):
        self.cm.extract_building_points(self.pc)
        self.assertIsInstance(self.cm._builder_cm, _pybuilder.CityModel)
        self.assertTrue(self.cm.extracted_points)
        b0 = self.cm._builder_cm.buildings[0]
        roof_pts = get_building_roof_points(b0)
        l, d = roof_pts.shape
        self.assertEqual(l, 220)
        self.assertEqual(d, 3)

    def test_ransac_filter(self):
        cm = CityModel(self.shp_file, p)
        cm.extract_building_points(self.pc)
        cm.building_points_RANSAC_outlier_remover(3, 250)
        self.assertIsInstance(cm._builder_cm, _pybuilder.CityModel)

    def test_statistical_outlier_filter(self):
        cm = CityModel(self.shp_file, p)
        cm.extract_building_points(self.pc, 1.0, 2.0)
        cm.building_points_statistical_outlier_remover(5, 1.5)
        self.assertIsInstance(cm._builder_cm, _pybuilder.CityModel)

    def test_compute_building_heights(self):
        cm = CityModel(self.shp_file, p)
        cm.extract_building_points(self.pc, 1.0, 2.0)
        dtm = ElevationModel(self.pc, 0.5, [2, 9])
        cm.compute_building_heights(dtm, 0.9, 0.95)
        self.assertIsInstance(self.cm._builder_cm, _pybuilder.CityModel)
        self.assertEqual(cm._builder_cm.buildings[0].height, 5)
        self.assertEqual(cm._builder_cm.buildings[3].height, 10)

    def test_readwrite_json(self):
        self.cm.to_JSON("tmpCM.json")
        cm2 = CityModel()
        cm2.from_JSON("tmpCM.json")
        self.assertIsInstance(cm2._builder_cm, _pybuilder.CityModel)
        os.unlink("tmpCM.json")

    def test_read_protobuf(self):
        pb_file = data_dir / "MinimalCase" / "PropertyMap.shp.pb"
        with open(pb_file, "rb") as src:
            pb_data = src.read()
        cm = CityModel()
        cm.load_protobuf(pb_data)
        self.assertIsInstance(cm._builder_cm, _pybuilder.CityModel)
        self.assertEqual(len(cm.buildings), 5)

    def test_convert_to_protobuf(self):
        pb_str = self.cm.to_protobuf()
        self.assertGreater(len(pb_str), 100)


if __name__ == "__main__":
    unittest.main()
