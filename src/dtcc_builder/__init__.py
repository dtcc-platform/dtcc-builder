from dtcc_builder import build
from dtcc_builder.Parameters import load_parameters
from dtcc_builder import builder_datamodel

# Add model extensions
from dtcc_model import CityModel, PointCloud

CityModel.add_processors(build.extract_buildingpoints, "extract_buildingpoints")
CityModel.add_processors(build.calculate_building_heights, "calculate_building_heights")


# Initialize logging
from dtcc_common import init_logging

init_logging("dtcc-builder")
