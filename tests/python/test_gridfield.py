import dtcc_io as io
from dtcc_model import Raster
import unittest
from dtcc_builder.model import raster_to_builder_gridfield
from pathlib import Path


data_dir = (Path(__file__).parent / "../data").resolve()


class TestGridField(unittest.TestCase):
    def test_gridfield_size(self):
        raster = io.load_raster(data_dir / "test_dem.tif")
        gridfield = raster_to_builder_gridfield(raster)
        self.assertEqual(gridfield.grid.xsize, 20)
        self.assertEqual(gridfield.grid.ysize, 40)
        self.assertEqual(gridfield.grid.xstep, 2.0)
        self.assertEqual(gridfield.grid.ystep, 2.0)

    def test_gridfield_values(self):
        raster = io.load_raster(data_dir / "test_dem.tif")
        gridfield = raster_to_builder_gridfield(raster)
        values = gridfield.values

        self.assertEqual(len(values), 20 * 40)

        # start in bottom left corner
        self.assertAlmostEqual(values[0], 76, places=5)
        self.assertAlmostEqual(values[1], 76.2, places=5)
        # top left corner
        self.assertAlmostEqual(values[(20 * 40) - 20], 0, places=5)
        # end in top right corner
        self.assertAlmostEqual(values[-1], 3.8, places=5)
