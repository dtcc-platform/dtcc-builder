from . import model
from . import builders
from . import parameters

from .builders import build, build_city, build_mesh, build_volume_mesh

# Add model extensions
from dtcc_model import City, PointCloud

City.add_processors(builders.compute_building_points, "compute_building_points")
City.add_processors(builders.compute_building_heights, "compute_building_heights")
