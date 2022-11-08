import unittest
import sys, os
from pathlib import Path

sys.path.append(str((Path(__file__).parent / "..").resolve()))

import _pybuilder
import Parameters
import PointCloud
import ElevationModel
import CityModel


data_dir = (Path(__file__).parent / "../../../unittests/data").resolve()
p = Parameters.load_parameters()


class TestCityModel(unittest.TestCase):
    def setUp(self):
        self.shp_file = data_dir / "MinimalCase" / "PropertyMap.shp"
        self.cm = CityModel.generate_citymodel(self.shp_file, p)
        self.pc = PointCloud.read_las_files(data_dir / "MinimalCase" / "pointcloud.las")

    def test_bounds(self):
        bbox = CityModel.building_bounds(self.shp_file)
        px, py, qx, qy = bbox
        self.assertAlmostEqual(px, -5.14247441879)
        self.assertAlmostEqual(qy, -1.09814696217)

    def test_padded_bounds(self):
        bbox = CityModel.building_bounds(self.shp_file, 10)
        px, py, qx, qy = bbox
        self.assertAlmostEqual(px, -15.14247441879)
        self.assertAlmostEqual(qy, 8.90185303782)

    def test_set_origin(self):
        bbox = CityModel.building_bounds(self.shp_file, 10)
        cm_with_origin = CityModel.set_origin(self.cm, bbox)
        self.assertIsInstance(cm_with_origin, _pybuilder.CityModel)

    def test_load(self):
        shp_file = data_dir / "MinimalCase" / "PropertyMap.shp"
        cm = CityModel.generate_citymodel(shp_file, p)
        self.assertIsInstance(cm, _pybuilder.CityModel)

    def test_clean(self):
        shp_file = data_dir / "MinimalCase" / "PropertyMap.shp"
        cm = CityModel.generate_citymodel(shp_file, p)
        clean_cm = CityModel.clean_citymodel(cm, 0.5)
        self.assertIsInstance(clean_cm, _pybuilder.CityModel)

    def test_extarct_points(self):
        cm_with_points = CityModel.extract_building_points(self.cm, self.pc, 1.0, 2.0)
        self.assertIsInstance(cm_with_points, _pybuilder.CityModel)
        b0 = cm_with_points.buildings[0]
        roof_pts = CityModel.get_building_roof_points(b0)
        l, d = roof_pts.shape
        self.assertEqual(l, 220)
        self.assertEqual(d, 3)

    def test_ransac_filter(self):
        cm_with_points = CityModel.extract_building_points(self.cm, self.pc, 1.0, 2.0)
        ransac_filtered = CityModel.building_points_RANSAC_outlier_remover(
            cm_with_points, 3, 250
        )
        self.assertIsInstance(ransac_filtered, _pybuilder.CityModel)

    def test_statistical_outlier_filter(self):
        cm_with_points = CityModel.extract_building_points(self.cm, self.pc, 1.0, 2.0)
        so_filtered = CityModel.building_points_statistical_outlier_remover(
            cm_with_points, 5, 1.5
        )
        self.assertIsInstance(so_filtered, _pybuilder.CityModel)

    def test_compute_building_heights(self):
        cm_with_points = CityModel.extract_building_points(self.cm, self.pc, 1.0, 2.0)
        dtm = ElevationModel.generate_elevation_model(self.pc, 0.5, [2, 9])
        cm_with_height = CityModel.compute_building_heights(
            cm_with_points, dtm, 0.9, 0.95
        )
        self.assertIsInstance(cm_with_height, _pybuilder.CityModel)
        self.assertEqual(cm_with_height.buildings[0].height, 5)
        self.assertEqual(cm_with_height.buildings[3].height, 10)

    def test_readwrite_json(self):
        CityModel.to_JSON(self.cm, "tmpCM.json")
        read_cm = CityModel.from_JSON("tmpCM.json")
        self.assertIsInstance(read_cm, _pybuilder.CityModel)
        os.unlink("tmpCM.json")

    def test_read_protobuf(self):
        pb_file = data_dir / "MinimalCase" / "PropertyMap.shp.pb"
        with open(pb_file, "rb") as src:
            pb_data = src.read()
        cm = CityModel.load_protobuf(pb_data)
        self.assertIsInstance(cm, _pybuilder.CityModel)
        self.assertEqual(len(cm.buildings), 5)

    def test_convert_to_protobuf(self):
        pb_str = CityModel.to_protobuf(self.cm)
        self.assertGreater(len(pb_str), 100)


if __name__ == "__main__":
    unittest.main()
