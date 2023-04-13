from dtcc_builder.scripts.dtcc_build import get_project_paths
import unittest
from pathlib import Path


class TestGetProjectPaths(unittest.TestCase):
    def test_get_project_paths_none(self):
        cwd = Path.cwd()
        self.assertEqual(get_project_paths(None), (cwd / "Parameters.json", cwd))

    def test_get_project_paths_dot(self):
        cwd = Path.cwd()
        self.assertEqual(get_project_paths("."), (cwd / "Parameters.json", cwd))

    def test_get_project_paths_file(self):
        wd = Path((Path(__file__).parent / "../data").resolve())
        self.assertEqual(get_project_paths(wd), (wd / "Parameters.json", wd))
