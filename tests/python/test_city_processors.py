import unittest
import dtcc_builder as builder
import dtcc_io as io
from pathlib import Path

project_dir = (Path(__file__).parent / "../data" / "MinimalCase").resolve()


class TestExtractPoints(unittest.TestCase):
    def test_extract_points(self):
        pc = io.load_pointcloud(project_dir / "pointcloud.las")
        city = io.load_city(project_dir / "PropertyMap.shp")

        city = city.extract_buildingpoints(pc)
        self.assertEqual(len(city.buildings[0].roofpoints), 216)
        self.assertEqual(len(city.buildings[3].roofpoints), 572)

    def test_calc_building_heights(self):
        pc = io.load_pointcloud(project_dir / "pointcloud.las")
        city = io.load_city(project_dir / "PropertyMap.shp")

        city = city.extract_buildingpoints(pc).calculate_building_heights()
        self.assertEqual(city.buildings[0].height, 5.0)
        self.assertEqual(city.buildings[3].height, 10.0)
