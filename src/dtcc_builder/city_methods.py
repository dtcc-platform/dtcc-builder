# Copyright (C) 2022 Dag Wästberg
# Licensed under the MIT License
#
# Modified by Anders Logg 2023
#
# This module provides functionality for city processing.

import numpy as np
from psutil import cpu_count
from math import ceil, log2
from time import time

import dtcc_model as model
from dtcc_model import City, PointCloud
from .logging import info, warning, error
from . import _dtcc_builder
from . import meshing
from . import model as builder_model
from . import parameters as builder_parameters


def compute_building_points(
    city: City,
    pointcloud: PointCloud,
    ground_margin=1.0,
    outlier_margin=2.0,
    statistical_outlier_remover=True,
    roof_outlier_neighbors=5,
    roof_outlier_margin=1.5,
    ransac_outlier_remover=False,
    ransac_outlier_margin=3.0,
    ransac_iterations=250,
):
    """
    Compute building points for the given city using point cloud data.

    Parameters
    ----------
    `city` : dtcc_model.City
        The city object for which to compute building points.
    `pointcloud` : dtcc_model.PointCloud
        The point cloud data associated with the city.
    `ground_margin` : float, optional
        The margin to use when filtering ground points, by default 1.0.
    `outlier_margin` : float, optional
        The margin to use when filtering outlier points, by default 2.0.
    `statistical_outlier_remover` : bool, optional
        Whether to use a statistical outlier remover, by default True.
    `roof_outlier_neighbors` : int, optional
        The number of neighbors to use when filtering roof outliers, by default 5.
    `roof_outlier_margin` : float, optional
        The margin to use when filtering roof outliers, by default 1.5.
    `ransac_outlier_remover` : bool, optional
        Whether to use a RANSAC outlier remover, by default False.
    `ransac_outlier_margin` : float, optional
        The margin to use when filtering RANSAC outliers, by default 3.0.
    `ransac_iterations` : int, optional
        The number of iterations to use when running RANSAC, by default 250.


    Returns
    -------
    `dtcc_model.City`
        The city object with computed building points.
    """
    info("Compute building points...")

    # Convert to builder model
    start_time = time()
    builder_city = builder_model.create_builder_city(city)
    builder_pointcloud = builder_model.create_builder_pointcloud(pointcloud)
    print(f"BBBB: Creating builder models took {time() - start_time} seconds")

    start_time = time()
    builder_pointcloud = _dtcc_builder.remove_vegetation(builder_pointcloud)
    print(f"BBBB: Removing vegetation took {time() - start_time} seconds")

    # Compute building points
    # builder_city = _dtcc_builder.compute_building_points(
    #     builder_city, builder_pointcloud, ground_margin, outlier_margin
    # )

    num_cores = cpu_count(logical=True)
    if num_cores is None or num_cores == 0:
        num_cores = 1
    num_tiles = ceil(log2(num_cores))
    num_tiles = max(2, num_tiles)
    num_tiles = min(8, num_tiles)
    info(f"Compute building points in parallel with {num_tiles}x{num_tiles} tiles")
    start_time = time()
    builder_city = _dtcc_builder.compute_building_points_parallel(
        builder_city,
        builder_pointcloud,
        ground_margin,
        outlier_margin,
        num_tiles,
        num_tiles,
    )
    print(f"BBBB: Computing building points took {time() - start_time} seconds")

    # Remove outliers
    start_time = time()
    if statistical_outlier_remover:
        builder_city = _dtcc_builder.remove_building_point_outliers_statistical(
            builder_city,
            roof_outlier_neighbors,
            roof_outlier_margin,
        )
    print(f"BBBB: Removing outliers took {time() - start_time} seconds")
    start_time = time()
    if ransac_outlier_remover:
        builder_city = _dtcc_builder.remove_building_point_outliers_ransac(
            builder_city,
            ransac_outlier_margin,
            ransac_iterations,
        )
    print(f"BBBB: Removing RANSAC outliers took {time() - start_time} seconds")
    # FIXME: Don't modify incoming data (city)

    # Convert back to city model
    start_time = time()
    roof_points = _dtcc_builder.building_roofpoints(builder_city)
    ground_points = _dtcc_builder.building_groundpoints(builder_city)
    print(f"BBBB: roof_points {len(roof_points)}")
    print(f"BBBB: roofpoint shapr {roof_points[0].shape}")
    for city_building, pts in zip(city.buildings, roof_points):
        city_building.roofpoints.points = pts
    for city_building, pts in zip(city.buildings, ground_points):
        ground_points = pts
        if len(ground_points) > 0:
            ground_z = ground_points[:, 2]
            city_building.ground_level = np.percentile(ground_z, 50)
    print(f"BBBB: Converting back to city model took {time() - start_time} seconds")
    return city


def compute_building_heights(
    city: model.City, min_building_height=2.5, roof_percentile=0.9
) -> model.City:
    """
    Compute building heights from roof points for the given city.

    Parameters
    ----------
    `city` : dtcc_model.City
        The city object for which to compute building heights.
    `parameters` : dict, optional
        A dictionary of parameters for the computation, by default None.

    Returns
    -------
    `dtcc_model.City`
        The city object with computed building heights.

    Notes
    -----
    - If a building's roof points are missing, the height is set to a minimum.
    - If a building's ground level height is missing, it is set to the ground height.
    - If a building's height is less than the minimum, it's set to the minimum height.


    """

    info("Computing building heights...")

    # Get parameters
    min_building_height = min_building_height
    roof_percentile = roof_percentile

    # FIXME: Don't modify incoming data (city)

    # Iterate over buildings
    num_too_low = 0
    for building in city.buildings:
        # Set building height to minimum height if points missing
        if len(building.roofpoints) == 0:
            # info(
            #     f"Building {building.uuid} has no roof points; setting height to minimum height f{min_building_height:.3f}m"
            # )
            building.height = min_building_height
            continue

        # Set ground level if missing
        if building.ground_level == 0 and len(city.terrain.shape) == 2:
            footprint_center = building.footprint.centroid
            ground_height = city.terrain.get_value(
                footprint_center.x, footprint_center.y
            )
            building.ground_level = ground_height

        # Calculate height
        z_values = building.roofpoints.points[:, 2]
        roof_top = np.percentile(z_values, roof_percentile * 100)
        height = roof_top - building.ground_level

        # Modify height if too small
        if height < min_building_height:
            height = min_building_height
            num_too_low += 1

        # Set building height
        building.height = height
    if num_too_low > 0:
        info(f"{num_too_low} buildings have height less than minimum height")

    return city


def extrude_buildings(city: City, mesh_resolution=5, zero_ground=False, cap_base=True):
    """
    Extrude buildings in the given city.

    Parameters
    ----------
    `city` : dtcc_model.City
        The city object for which to extrude buildings.
    `zero_ground` : bool, optional
        Whether to set the ground level to zero, by default False.
    `cap_base` : bool, optional
        Whether to cap the base of the buildings, by default True.

    Returns
    -------
    `dtcc_model.City`
        The city object with building.mesh set to the extrusion.
    """
    info("Extruding buildings...")
    for b in city.buildings:
        extrusion = meshing.extrude_building(b, mesh_resolution, zero_ground, cap_base)
        b.mesh = extrusion
    return city
