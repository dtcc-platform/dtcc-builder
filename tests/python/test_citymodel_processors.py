import unittest
import dtcc_builder as builder
import dtcc_io as io
from pathlib import Path

project_dir = (Path(__file__).parent / "../data" / "MinimalCase").resolve()


class TestExtractPoints(unittest.TestCase):
    def test_extract_points(self):
        pc = io.load_pointcloud(project_dir / "pointcloud.las")
        cm = io.load_citymodel(project_dir / "PropertyMap.shp")

        cm = cm.extract_buildingpoints(pc)
        self.assertEqual(len(cm.buildings[0].roofpoints), 216)
        self.assertEqual(len(cm.buildings[3].roofpoints), 572)

    def test_calc_building_heights(self):
        pc = io.load_pointcloud(project_dir / "pointcloud.las")
        cm = io.load_citymodel(project_dir / "PropertyMap.shp")

        cm = cm.extract_buildingpoints(pc).calculate_building_heights()
        self.assertEqual(cm.buildings[0].height, 5.0)
        self.assertEqual(cm.buildings[3].height, 10.0)
