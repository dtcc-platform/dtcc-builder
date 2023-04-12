import unittest
from pathlib import Path
import sys
from dtcc_model import CityModel, GridField2D, Mesh, VolumeMesh

from dtcc_builder import Build, load_parameters

data_dir = (Path(__file__).parent / "../data").resolve()

footprint_path = data_dir / "MinimalCase" / "PropertyMap.shp"
pointcloud_path = data_dir / "MinimalCase" / "pointcloud.las"
print(footprint_path.is_file())


class TestBuild(unittest.TestCase):
    def test_build_citymodel(self):
        p = load_parameters()
        print(footprint_path)
        cm, dtm = Build.build_citymodel(footprint_path, pointcloud_path, p)
        self.assertIsInstance(cm, CityModel)
        # self.assertIsInstance(dtm, ElevationModel)

    def test_build_mesh(self):
        p = load_parameters()
        p["Debug"] = False
        volume_mesh, surface_mesh = Build.build_volume_mesh(
            footprint_path, pointcloud_path, p
        )
        self.assertIsInstance(volume_mesh, VolumeMesh)
        self.assertIsInstance(surface_mesh, Mesh)
