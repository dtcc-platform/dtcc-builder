from dtcc_builder import build
from dtcc_builder.Parameters import load_parameters
from dtcc_builder import builder_datamodel

# from dtcc_builder.Build import build_citymodel, build_surface_meshes, build_volume_mesh

from dtcc_model import CityModel, PointCloud

CityModel.add_processors(build.extract_buildingpoints, "extract_buildingpoints")
CityModel.add_processors(build.calculate_building_heights, "calculate_building_heights")
