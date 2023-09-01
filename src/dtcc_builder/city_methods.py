# Copyright (C) 2022 Dag Wästberg
# Licensed under the MIT License
#
# Modified by Anders Logg 2023
#
# This module provides functionality for city processing.

import numpy as np

import dtcc_model as model
from .logging import info, warning, error
from . import _dtcc_builder
from . import model as builder_model
from . import parameters as builder_parameters


def compute_building_points(city, pointcloud, parameters: dict = None):
    """
    Compute building points for the given city using point cloud data.

    Parameters
    ----------
    `city` : dtcc_model.City
        The city object for which to compute building points.
    `pointcloud` : dtcc_model.PointCloud
        The point cloud data associated with the city.
    `parameters` : dict, optional
        A dictionary of parameters for the computation, by default None.

    Returns
    -------
    `dtcc_model.City`
        The city object with computed building points.
    """
    info("Compute building points...")

    # Get parameters
    p = parameters or builder_parameters.default()

    pointcloud = _dtcc_builder.pointcloud.remove_vegetation()
    # Convert to builder model
    builder_city = builder_model.create_builder_city(city)
    builder_pointcloud = builder_model.create_builder_pointcloud(pointcloud)

    # Compute building points
    builder_city = _dtcc_builder.compute_building_points(
        builder_city, builder_pointcloud, p["ground_margin"], p["outlier_margin"]
    )

    # Remove outliers
    if p["statistical_outlier_remover"]:
        builder_city = _dtcc_builder.remove_building_point_outliers_statistical(
            builder_city,
            p["roof_outlier_neighbors"],
            p["roof_outlier_margin"],
        )
    if p["ransac_outlier_remover"]:
        builder_city = _dtcc_builder.remove_building_point_outliers_ransac(
            builder_city,
            p["ransac_outlier_margin"],
            p["ransac_iterations"],
        )

    # FIXME: Don't modify incoming data (city)

    # Convert back to city model
    for city_building, builder_buildings in zip(city.buildings, builder_city.buildings):
        city_building.roofpoints.points = np.array(
            [[p.x, p.y, p.z] for p in builder_buildings.roof_points]
        )
        ground_points = np.array(
            [[p.x, p.y, p.z] for p in builder_buildings.ground_points]
        )
        if len(ground_points) > 0:
            ground_z = ground_points[:, 2]
            city_building.ground_level = np.percentile(ground_z, 50)

    return city


def compute_building_heights(city: model.City, parameters: dict = None) -> model.City:
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
    p = parameters or builder_parameters.default()
    min_building_height = p["min_building_height"]
    roof_percentile = p["roof_percentile"]

    # FIXME: Don't modify incoming data (city)

    # Iterate over buildings
    for building in city.buildings:
        # Set building height to minimum height if points missing
        if len(building.roofpoints) == 0:
            info(
                f"Building {building.uuid} has no roof points; setting height to minimum height f{min_building_height:.3f}m"
            )
            building.height = min_building_height
            continue

        # Set ground level if missing
        if building.ground_level == 0:
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
            info(
                f"Building {building.uuid} to low ({height:.3f}m); setting height to minimum height f{min_building_height:.3f}m"
            )
            height = min_building_height

        # Set building height
        building.height = height

    return city
