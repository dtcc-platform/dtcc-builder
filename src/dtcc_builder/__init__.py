from . import model
from . import builders
from . import parameters
from . import city_methods

from .builders import (
    compute_domain_bounds,
    build,
    build_city,
    build_mesh,
    build_volume_mesh,
)

# Add model extensions
from dtcc_model import City, PointCloud

City.add_methods(city_methods.compute_building_points, "compute_building_points")
City.add_methods(city_methods.compute_building_heights, "compute_building_heights")

__all__ = [
    "compute_domain_bounds",
    "build",
    "build_city",
    "build_mesh",
    "build_volume_mesh",
]
