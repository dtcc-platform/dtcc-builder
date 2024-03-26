from dtcc_builder.geometry_builders.surface import extrude_surface
import unittest
import numpy as np

from dtcc_model import Building, MultiSurface, Surface


class TestExtrudeSurface(unittest.TestCase):
    def test_extrude_simple_surface(self):
        surface = Surface(
            vertices=np.array([[0, 0, 10], [10, 0, 10], [10, 10, 10], [0, 10, 10]])
        )

        extruded = extrude_surface(surface, 2)
        self.assertEqual(len(extruded.surfaces), 6)
        self.assertAlmostEqual(extruded.surfaces[2].vertices[:, 2].max(), 10)
        self.assertAlmostEqual(extruded.surfaces[2].vertices[:, 2].min(), 2)


if __name__ == "__main__":
    unittest.main()
