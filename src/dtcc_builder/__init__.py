from . import methods
from . import model
from . import parameters

from .methods import build

# Add model extensions
from dtcc_model import City, PointCloud

City.add_processors(methods.extract_buildingpoints, "extract_buildingpoints")
City.add_processors(methods.calculate_building_heights, "calculate_building_heights")


# Initialize logging
from dtcc_common import init_logging

init_logging("dtcc-builder")
