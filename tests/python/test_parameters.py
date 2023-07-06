import dtcc_builder

import unittest


class TestParameters(unittest.TestCase):
    def test_default_parameters(self):
        p = dtcc_builder.parameters.default()
        self.assertIsInstance(p, dict)
        self.assertEqual(p["model_name"], "DTCC")
