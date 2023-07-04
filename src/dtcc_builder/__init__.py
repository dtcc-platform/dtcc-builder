from . import model
from . import builders
from . import parameters
from . import city_processing

from .builders import build, build_dem, build_city, build_mesh, build_volume_mesh

# Add model extensions
from dtcc_model import City, PointCloud

City.add_processors(city_processing.compute_building_points, "compute_building_points")
City.add_processors(
    city_processing.compute_building_heights, "compute_building_heights"
)
