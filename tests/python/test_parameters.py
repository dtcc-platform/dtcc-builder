from dtcc_builder.Parameters import load_parameters

from pathlib import Path
import unittest

project_dir = (Path(__file__).parent / "../data" / "MinimalCase").resolve()

print(project_dir)


class TestLoadParameters(unittest.TestCase):
    def test_load_parameters(self):
        p = load_parameters()
        self.assertIsInstance(p, dict)
        self.assertEqual(p["model-name"], "DTCC")

    def test_load_parameters_file(self):
        p = load_parameters(project_dir / "Parameters.json")
        self.assertIsInstance(p, dict)
        print(p)
        self.assertEqual(p["model-name"], "MinimalCase")
