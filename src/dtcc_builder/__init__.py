from . import build
from . import builder_datamodel
from . import parameters

# Add model extensions
from dtcc_model import City, PointCloud

City.add_processors(build.extract_buildingpoints, "extract_buildingpoints")
City.add_processors(build.calculate_building_heights, "calculate_building_heights")


# Initialize logging
from dtcc_common import init_logging

init_logging("dtcc-builder")
