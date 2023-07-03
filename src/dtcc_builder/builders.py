# Copyright (C) 2022 Dag Wästberg
# Licensed under the MIT License
#
# Modified by Anders Logg 2023
#
# This module provides the main methods of DTCC Builder.

import numpy as np
import rasterio
import pypoints2grid
from pathlib import Path
from typing import Tuple, List

import dtcc_model as model
import dtcc_io as io
from .logging import info, warning, error
from . import _dtcc_builder
from . import model as builder_model


def compute_domain_bounds(params_or_footprint, las_path=None, params=None):
    "Compute domain bounds from footprint and pointcloud"

    info("Computing domain bounds...")

    if las_path is None and params is None:
        p = params_or_footprint
        building_footprint_path = p["data-directory"] / p["buildings-filename"]
        pointcloud_path = p["pointcloud-directory"]
    elif (
        params_or_footprint is not None and las_path is not None and params is not None
    ):
        p = params
        building_footprint_path = params_or_footprint
        pointcloud_path = las_path
    else:
        raise ValueError(
            "must provide either params only or footprint and las_path and params"
        )
    if p["auto_domain"]:
        footprint_bounds = io.city.building_bounds(
            building_footprint_path, p["domain_margin"]
        )
        las_bounds = io.pointcloud.calc_las_bounds(pointcloud_path)
        domain_bounds = io.bounds.bounds_intersect(footprint_bounds, las_bounds)
        origin = domain_bounds[:2]
        p["x0"] = origin[0]
        p["y0"] = origin[1]
        p["x_min"] = 0.0
        p["y_min"] = 0.0
        p["x_max"] = domain_bounds[2] - domain_bounds[0]
        p["y_max"] = domain_bounds[3] - domain_bounds[1]
    else:
        origin = (p["x0"], p["y0"])
        domain_bounds = (
            p["x0"] + p["x_min"],
            p["y0"] + p["y_min"],
            p["x0"] + p["x_max"],
            p["y0"] + p["y_max"],
        )

    domain_bounds = model.Bounds(
        xmin=domain_bounds[0],
        ymin=domain_bounds[1],
        xmax=domain_bounds[2],
        ymax=domain_bounds[3],
    )
    return (origin, domain_bounds)


def build_dem(
    pointcloud: model.PointCloud, bounds, cell_size: float, window_size: int = 3
) -> model.Raster:
    info("Building DEM...")

    if (
        len(pointcloud.classification) == len(pointcloud.points)
    ) and 2 in pointcloud.used_classifications():
        ground_point_idx = np.where(np.isin(pointcloud.classification, [2, 9]))[0]
        ground_points = pointcloud.points[ground_point_idx]
    else:
        ground_points = pointcloud.points
    print(f"generating dem with bounds {bounds.tuple}")
    dem = pypoints2grid.points2grid(
        ground_points, cell_size, bounds.tuple, window_size=window_size
    )

    print("dem shape: ", dem.shape)
    print("cell size: ", cell_size)
    dem_raster = model.Raster()
    dem_raster.data = dem
    dem_raster.nodata = 0
    dem_raster.georef = rasterio.transform.from_origin(
        bounds.west, bounds.north, cell_size, cell_size
    )

    print(f"generated dem with bounds {dem_raster.bounds}")
    print(f"generated dem with georef \n{dem_raster.georef}")
    dem_raster = dem_raster.fill_holes()

    return dem_raster


def compute_building_heights(
    city, roof_percentile=0.9, min_height=2.5, overwrite=False
):
    info("Computing building heights...")

    for building in city.buildings:
        if building.height > 0 and not overwrite:
            continue
        if len(building.roofpoints) == 0:
            building.height = min_height
            continue
        z_values = building.roofpoints.points[:, 2]
        roof_top = np.percentile(z_values, roof_percentile * 100)

        if building.ground_level == 0:
            if len(city.terrain.shape) == 2:
                footprint_center = building.footprint.centroid
                ground_height = city.terrain.get_value(
                    footprint_center.x, footprint_center.y
                )
                building.ground_level = ground_height
        height = roof_top - building.ground_level
        if height < min_height:
            height = min_height
        building.height = height

    return city


def compute_building_points(
    city,
    pointcloud,
    ground_padding=2.0,
    ground_outlier_margin=1,
    roof_outlier_margin=1.5,
    roof_outlier_neighbors=5,
    roof_ransac_margin=3.0,
    roof_ransac_iterations=150,
):
    info("Compute building points...")

    builder_city = builder_model.create_builder_city(city)
    builder_pointcloud = builder_model.create_builder_pointcloud(pointcloud)
    builder_city = _dtcc_builder.extractRoofPoints(
        builder_city,
        builder_pointcloud,
        ground_padding,
        ground_outlier_margin,
        roof_outlier_margin,
        roof_outlier_neighbors,
        roof_ransac_margin,
        roof_ransac_iterations,
    )
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


def build_city(
    city: model.City,
    point_cloud: model.PointCloud,
    bounds: model.Bounds,
    parameters: dict,
) -> model.City:
    """
    Build city from building footprints.

    Developer note: Consider introducing a new class named Footprints
    so that a city can be built from footprints and point cloud data.
    It is somewhat strange that the input to this function is a city.
    """

    info("Building city...")

    # Shortcut
    p = parameters

    # Remove outliers from point cloud
    point_cloud = point_cloud.remove_global_outliers(p["outlier_margin"])

    # Build elevation model FIXME: Don't modify the input city
    city.terrain = build_dem(point_cloud, bounds, p["elevation_model_resolution"])

    # Compute building points
    city = compute_building_points(
        city,
        point_cloud,
        p["ground_margin"],
        p["outlier_margin"],
        p["roof_outlier_margin"],
        p["outlier_neighbors"],
        p["ransac_outlier_margin"] if p["ransac_outlier_remover"] else 0,
        p["ransac_iterations"] if p["ransac_outlier_remover"] else 0,
    )

    # Compute building heights
    city = compute_building_heights(
        city, p["roof_percentile"], p["min_building_height"], overwrite=True
    )

    return city


def build_mesh(
    city: model.City, parameters: dict
) -> Tuple[model.Mesh, List[model.Mesh]]:
    """
    Build mesh for city.

    This function builds a mesh for the city, returning two different
    meshes: one for the ground and one for the buildings.
    """

    info("Building meshes for city...")

    # Shortcut
    p = parameters

    # Convert to builder model
    builder_city = builder_model.create_builder_city(city)
    builder_dem = builder_model.raster_to_builder_gridfield(city.terrain)
    bounds = (
        city.bounds.xmin,
        city.bounds.ymin,
        city.bounds.xmax,
        city.bounds.ymax,
    )

    # Simplify city
    simple_city = _dtcc_builder.SimplifyCity(
        builder_city, bounds, p["min_building_distance"], p["min_vertex_distance"]
    )

    # Build meshes
    meshes = _dtcc_builder.BuildSurface3D(
        simple_city, builder_dem, p["mesh_resolution"]
    )

    # Extract meshes and merge building meshes
    ground_mesh = meshes[0]
    building_meshes = meshes[1:]
    building_meshes = _dtcc_builder.MergeSurfaces3D(building_meshes)

    # Convert back to DTCC model
    dtcc_ground_mesh = builder_model.builder_mesh_to_mesh(ground_mesh)
    dtcc_building_mesh = builder_model.builder_mesh_to_mesh(building_meshes)

    return dtcc_ground_mesh, dtcc_building_mesh


def _debug(mesh, step, p):
    "Debug volume meshing"
    if not p["debug"]:
        return
    if isinstance(mesh, _dtcc_builder.Mesh):
        mesh = builder_model.builder_mesh_to_mesh(mesh)
    else:
        mesh = builder_model.builder_volume_mesh_to_volume_mesh(mesh)
    mesh.save(p["output_directory"] / f"mesh_step{step}.vtu")


def build_volume_mesh(
    city: model.City, parameters: dict
) -> Tuple[model.VolumeMesh, model.Mesh]:
    """Build volume mesh for city.

    This function builds a boundary conforming volume mesh for the city,
    returning both the volume mesh and its corresponding  boundary mesh.
    """

    info("Building volume mesh for city...")

    # Shortcut
    p = parameters

    # Convert to builder model
    builder_city = builder_model.create_builder_city(city)
    builder_dem = builder_model.raster_to_builder_gridfield(city.terrain)
    bounds = (
        city.bounds.xmin,
        city.bounds.ymin,
        city.bounds.xmax,
        city.bounds.ymax,
    )

    # Simplify city
    simple_city = _dtcc_builder.SimplifyCity(
        builder_city, bounds, p["min_building_distance"], p["min_vertex_distance"]
    )

    # Step 3.1: Build 2D mesh
    mesh = _dtcc_builder.BuildMesh2D(
        simple_city,
        bounds,
        p["mesh_resolution"],
    )
    _debug(mesh, "3.1", p)

    # Step 3.2: Build 3D mesh (layer 3D mesh)
    volume_mesh = _dtcc_builder.BuildVolumeMesh(
        mesh, p["domain_height"], p["mesh_resolution"]
    )
    _debug(volume_mesh, "3.2", p)

    # Step 3.3: Smooth 3D mesh (set ground height)
    top_height = p["domain_height"] + city.terrain.data.mean()
    volume_mesh = _dtcc_builder.smooth_volume_mesh(
        volume_mesh,
        simple_city,
        builder_dem,
        top_height,
        False,
        p["smoothing_max_iterations"],
        p["smoothing_relative_tolerance"],
    )
    _debug(volume_mesh, "3.3", p)

    # Step 3.4: Trim 3D mesh (remove building interiors)
    volume_mesh = _dtcc_builder.TrimVolumeMesh(volume_mesh, mesh, simple_city)
    _debug(volume_mesh, "3.4", p)

    # Step 3.5: Smooth 3D mesh (set ground and building heights)
    volume_mesh = _dtcc_builder.smooth_volume_mesh(
        volume_mesh,
        simple_city,
        builder_dem,
        top_height,
        True,
        p["smoothing_max_iterations"],
        p["smoothing_relative_tolerance"],
    )
    _debug(volume_mesh, "3.5", p)

    # Compute boundary mesh
    volume_mesh_boundary = _dtcc_builder.ExtractBoundary3D(volume_mesh)

    # Convert back to DTCC model
    dtcc_volume_mesh = builder_model.builder_volume_mesh_to_volume_mesh(volume_mesh)
    dtcc_volume_mesh_boundary = builder_model.builder_mesh_to_mesh(volume_mesh_boundary)

    return dtcc_volume_mesh, dtcc_volume_mesh_boundary


def build(parameters):
    """
    Build city and city meshes.

    This function reads data from the specified data directory
    and builds a city and its corresponding meshes. The same
    thing can be accomplished by calling the individual build_*
    functions, but this function is provided as a convenience
    and takes care of loading and saving data to files.
    """

    # Shortcut
    p = parameters

    # Get paths for input data
    buildings_path = p["data_directory"] / p["buildings_filename"]
    pointcloud_path = p["pointcloud_directory"]
    if not buildings_path.exists():
        error(f"Unable to read buildings file {buildings_path}")
    if not pointcloud_path.exists():
        error(f"Unable to read pointcloud directory {pointcloud_path}")

    # Compute domain bounds
    origin, bounds = compute_domain_bounds(buildings_path, pointcloud_path, p)
    info(bounds)

    # Load city
    city = io.load_city(
        buildings_path,
        uuid_field=p["uuid_field"],
        height_field=p["height_field"],
        bounds=bounds,
    )

    # Load point cloud
    point_cloud = io.load_pointcloud(pointcloud_path, bounds=bounds)

    # Build city
    build_city(city, point_cloud, bounds, p)

    # Save city to file
    city_name = p["output_directory"] / "city"
    if p["save_protobuf"]:
        io.save_city(city, city_name.with_suffix(".pb"))
    if p["save_shp"]:
        io.save_city(city, city_name.with_suffix(".shp"))
    if p["save_json"]:
        io.save_city(city, city_name.with_suffix(".json"))

    # Build mesh
    if p["build_mesh"]:
        ground_mesh, building_mesh = build_mesh(city, p)

        # Save meshes to file
        ground_mesh_name = p["output_directory"] / "ground_mesh"
        building_mesh_name = p["output_directory"] / "building_mesh"
        if p["save_protobuf"]:
            ground_mesh.save(ground_mesh_name.with_suffix(".pb"))
            building_mesh.save(building_mesh_name.with_suffix(".pb"))
        if p["save_vtk"]:
            ground_mesh.save(ground_mesh_name.with_suffix(".vtk"))
            building_mesh.save(building_mesh_name.with_suffix(".vtk"))
        if p["save_stl"]:
            ground_mesh.save(ground_mesh_name.with_suffix(".stl"))
            building_mesh.save(building_mesh_name.with_suffix(".stl"))
        if p["save_obj"]:
            ground_mesh.save(ground_mesh_name.with_suffix(".obj"))
            building_mesh.save(building_mesh_name.with_suffix(".obj"))

    # Build volume mesh
    if p["build_volume_mesh"]:
        mesh, volume_mesh = build_volume_mesh(city, p)

        # Save meshes to file
        mesh_name = p["output_directory"] / "mesh"
        volume_mesh_name = p["output_directory"] / "volume_mesh"
        if p["save_protobuf"]:
            mesh.save(mesh_name.with_suffix(".pb"))
            volume_mesh.save(volume_mesh_name.with_suffix(".pb"))
        if p["save_vtk"]:
            mesh.save(mesh_name.with_suffix(".vtk"))
            volume_mesh.save(volume_mesh_name.with_suffix(".vtk"))
        if p["save_stl"]:
            mesh.save(mesh_name.with_suffix(".stl"))
        if p["save_obj"]:
            mesh.save(mesh_name.with_suffix(".obj"))
