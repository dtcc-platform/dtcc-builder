from . import model
from . import builders
from . import parameters
from . import city_methods
from . import meshing

from .builders import (
    build,
    build_city,
    calculate_bounds,
    build_terrain_mesh,
    build_building_meshes,
    build_city_surface_mesh,
    build_volume_mesh,
)

# Add model extensions
from dtcc_model import City, PointCloud

City.add_methods(city_methods.compute_building_points, "compute_building_points")
City.add_methods(city_methods.compute_building_heights, "compute_building_heights")

__all__ = [
    "build",
    "build_city",
    "build_terrain_mesh",
    "build_building_meshes",
    "build_volume_mesh",
    "build_city_surface_mesh",
    "calculate_bounds",
]
