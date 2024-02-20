from dtcc_model import PointCloud, Terrain, Raster, Mesh, Surface

from dtcc_builder.model import (
    raster_to_builder_gridfield,
    builder_mesh_to_mesh,
    create_builder_polygon,
)

import numpy as np
from pypoints2grid import points2grid
from affine import Affine


def build_terrain_mesh(
    dem: Raster = None,
    pointcloud: PointCloud = None,
    subdomains: [Surface] = None,
    max_mesh_size=10,
    min_mesh_angle=25,
    smoothing=3,
) -> Mesh:
    if dem is None and pointcloud is None:
        raise ValueError("Either dem raster or pointcloud must be provided.")
    if dem is None:
        dem = build_terrain_raster(pointcloud, cell_size=max_mesh_size / 2)
    _builder_gridfield = raster_to_builder_gridfield(dem)
    if subdomains is None:
        subdomains = []
    else:
        subdomains = [create_builder_polygon(sub.to_polygon) for sub in subdomains]
    terrain_mesh = _dtcc_builder.build_terrain_mesh(
        _builder_gridfield, max_mesh_size, min_mesh_angle, smoothing, subdomains
    )


def build_terrain_raster(
    pc: PointCloud, cell_size, bounds=None, window_size=3, radius=0, ground_only=True
) -> Raster:
    """
    Rasterize a point cloud into a `Raster` object.

    Args:
        cell_size (float): The size of the raster cells in meters.
        bounds (Bounds): The bounds of the area to rasterize (default None, uses the bounds of the point cloud).
        window_size (int): The size of the window for the interpolation (default 3).
        radius (float): The radius of the search for the interpolation (default 0).
        ground_only (bool): Whether to only use ground points for the rasterization (default True).

    Returns:
        Raster: A `Raster` object representing the rasterized point cloud.
    """
    if (
        ground_only
        and (len(pc.classification) == len(pc.points))
        and 2 in pc.used_classifications()
    ):
        ground_point_idx = np.where(np.isin(pc.classification, [2, 9]))[0]
        ground_points = pc.points[ground_point_idx]
    else:
        ground_points = pc.points
    if bounds is None:
        if pc.bounds is None or pc.bounds.area == 0:
            pc.calculate_bounds()
        bounds = pc.bounds
    dem = points2grid(
        ground_points, cell_size, bounds.tuple, window_size=window_size, radius=radius
    )
    dem_raster = Raster()
    dem_raster.data = dem
    dem_raster.nodata = 0
    dem_raster.georef = Affine.translation(bounds.minx, bounds.maxy) * Affine.scale(
        cell_size, -cell_size
    )
    dem_raster = dem_raster.fill_holes()
    return dem_raster
