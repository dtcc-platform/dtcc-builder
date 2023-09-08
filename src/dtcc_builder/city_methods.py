# Copyright (C) 2022 Dag Wästberg
# Licensed under the MIT License
#
# Modified by Anders Logg 2023
#
# This module provides functionality for city processing.

import numpy as np

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
    builder_city = builder_model.create_builder_city(city)
    builder_pointcloud = builder_model.create_builder_pointcloud(pointcloud)

    builder_pointcloud = _dtcc_builder.remove_vegetation(builder_pointcloud)

    # Compute building points
    builder_city = _dtcc_builder.compute_building_points(
        builder_city, builder_pointcloud, ground_margin, outlier_margin
    )

    # Remove outliers
    if statistical_outlier_remover:
        builder_city = _dtcc_builder.remove_building_point_outliers_statistical(
            builder_city,
            roof_outlier_neighbors,
            roof_outlier_margin,
        )
    if ransac_outlier_remover:
        builder_city = _dtcc_builder.remove_building_point_outliers_ransac(
            builder_city,
            ransac_outlier_margin,
            ransac_iterations,
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
    for building in city.buildings:
        # Set building height to minimum height if points missing
        if len(building.roofpoints) == 0:
            info(
                f"Building {building.uuid} has no roof points; setting height to minimum height f{min_building_height:.3f}m"
            )
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
            info(
                f"Building {building.uuid} to low ({height:.3f}m); setting height to minimum height f{min_building_height:.3f}m"
            )
            height = min_building_height

        # Set building height
        building.height = height

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
