import unittest
import sys
from pathlib import Path
sys.path.append( str((Path(__file__).parent / "..").resolve() ))

import _pybuilder
import Parameters
import CityModelGenerator

data_dir = (Path(__file__).parent / "../../unittests/data").resolve()
p = Parameters.loadParameters()

class TestCityModel(unittest.TestCase):

    def test_load(self):
        shp_file = data_dir / "MinimalCase" / "PropertyMap.shp"
        cm = CityModelGenerator.GenerateCityModel(shp_file,p)
        self.assertIsInstance(cm, _pybuilder.CityModel)

if __name__ == '__main__':
    unittest.main()