# Copyright (C) 2022 Dag Wästberg
# Licensed under the MIT License
#
# Modified by Anders Logg 2023
#
# This module provides the main methods of DTCC Builder.

import numpy as np
import rasterio.transform
import pypoints2grid
from pathlib import Path
from typing import Tuple

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


def build_mesh(
    city: model.City,
    mesh_resolution,
    domain_height,
    min_building_dist,
    min_vertex_dist,
    debug=False,
) -> Tuple[model.Mesh, model.VolumeMesh]:
    """Build boundary conforming mesh and volume mesh for a city"""

    info("Building boundary conforming city mesh...")

    builder_city = builder_model.create_builder_city(city)
    bounds = (
        city.bounds.xmin,
        city.bounds.ymin,
        city.bounds.xmax,
        city.bounds.ymax,
    )
    print(f"Building mesh with bounds {bounds}")
    print(f"dem bounds {city.terrain.bounds}")
    simple_city = _dtcc_builder.SimplifyCity(
        builder_city, bounds, min_building_dist, min_vertex_dist
    )
    # simple_city = _dtcc_builder.CleanCity(simple_city, min_vertex_dist)

    builder_dem = builder_model.raster_to_builder_gridfield(city.terrain)

    # Step 3.1: Build 2D mesh
    mesh_2D = _dtcc_builder.BuildMesh2D(
        simple_city,
        bounds,
        mesh_resolution,
    )

    if debug:
        builder_model.builder_mesh_to_mesh(mesh_2D).save("mesh_step3.1.vtu")

    # Step 3.2: Build 3D mesh (layer 3D mesh)
    mesh_3D = _dtcc_builder.BuildVolumeMesh(mesh_2D, domain_height, mesh_resolution)
    if debug:
        builder_model.builder_volume_mesh_to_volume_mesh(mesh_3D).save(
            "meshing_step3.2.vtu"
        )

    # FIXME: Make parameters
    max_iterations = 1000
    relative_tolerance = 1e-3

    # Step 3.3: Smooth 3D mesh (set ground height)
    top_height = domain_height + city.terrain.data.mean()
    mesh_3D = _dtcc_builder.smooth_volume_mesh(
        mesh_3D,
        simple_city,
        builder_dem,
        top_height,
        False,
        max_iterations,
        relative_tolerance,
    )

    if debug:
        builder_model.builder_volume_mesh_to_volume_mesh(mesh_3D).save(
            "meshing_step3.3.vtu"
        )

    # Step 3.4: Trim 3D mesh (remove building interiors)
    mesh_3D = _dtcc_builder.TrimVolumeMesh(mesh_3D, mesh_2D, simple_city)
    if debug:
        builder_model.builder_volume_mesh_to_volume_mesh(mesh_3D).save(
            "meshing_step3.4.vtu"
        )

    # Step 3.5: Smooth 3D mesh (set ground and building heights)
    mesh_3D = _dtcc_builder.smooth_volume_mesh(
        mesh_3D,
        simple_city,
        builder_dem,
        top_height,
        True,
        max_iterations,
        relative_tolerance,
    )

    if debug:
        builder_datamodel.builder_volume_mesh_to_volume_mesh(mesh_3D).save(
            "meshing_step3.5.vtu"
        )
    surface_3d = _dtcc_builder.ExtractBoundary3D(mesh_3D)

    # Convert back to DTCC model
    dtcc_mesh = builder_model.builder_mesh_to_mesh(surface_3d)
    dtcc_volume_mesh = builder_model.builder_volume_mesh_to_volume_mesh(mesh_3D)

    return dtcc_mesh, dtcc_volume_mesh


def build_meshes(city: model.City, min_building_dist, min_vertex_dist, mesh_resolution):
    """Build non boundary conforming meshes for a city"""

    info("Building city meshes (non boundary conforming)...")

    bounds = (
        city.bounds.xmin,
        city.bounds.ymin,
        city.bounds.xmax,
        city.bounds.ymax,
    )
    builder_city = builder_model.create_builder_city(city)
    simple_city = _dtcc_builder.SimplifyCity(
        builder_city, bounds, min_building_dist, min_vertex_dist
    )
    # simple_city = _dtcc_builder.CleanCity(simple_city, min_vertex_dist)

    builder_dem = builder_model.raster_to_builder_gridfield(city.terrain)

    surfaces = _dtcc_builder.BuildSurface3D(simple_city, builder_dem, mesh_resolution)

    ground_surface = surfaces[0]
    building_surfaces = surfaces[1:]
    building_surfaces = _dtcc_builder.MergeSurfaces3D(building_surfaces)

    dtcc_ground_surface = builder_model.builder_mesh_to_mesh(ground_surface)
    dtcc_building_surfaces = builder_model.builder_mesh_to_mesh(building_surfaces)

    return dtcc_ground_surface, dtcc_building_surfaces


def build(parameters, city_only=False, mesh_only=False):
    """Build city and/or volume mesh according to given parameters"""

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

    # Load point cloud and remove outliers
    point_cloud = io.load_pointcloud(pointcloud_path, bounds=bounds)
    point_cloud = point_cloud.remove_global_outliers(p["outlier_margin"])

    # Build elevation model
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

    # Save city to file
    prefix = "city"
    if p["save_protobuf"]:
        io.save_city(city, p["output_directory"] / f"{prefix}.pb")
    if p["save_shp"]:
        io.save_city(city, p["output_directory"] / f"{prefix}.shp")
    if p["save_json"]:
        io.save_city(city, p["output_directory"] / f"{prefix}.json")

    # Build conforming mesh
    mesh, volume_mesh = build_mesh(
        city,
        p["mesh_resolution"],
        p["domain_height"],
        p["min_building_distance"],
        p["min_vertex_distance"],
        p["debug"],
    )

    # Build non-conforming meshes
    ground_mesh, building_mesh = build_meshes(
        city,
        p["min_building_distance"],
        p["min_vertex_distance"],
        p["mesh_resolution"],
    )

    # Save meshes to file
    prefix = "mesh"
    if p["save_protobuf"]:
        mesh.save(p["output_directory"] / f"city_{prefix}.pb")
        volume_mesh.save(p["output_directory"] / f"city_volume_{prefix}.pb")
        ground_mesh.save(p["output_directory"] / f"ground_{prefix}.pb")
        building_mesh.save(p["output_directory"] / f"building_{prefix}.pb")
    if p["save_vtk"]:
        mesh.save(p["output_directory"] / f"city_{prefix}.vtk")
        volume_mesh.save(p["output_directory"] / f"city_volume_{prefix}.vtk")
        ground_mesh.save(p["output_directory"] / f"ground_{prefix}.vtk")
        building_mesh.save(p["output_directory"] / f"building_{prefix}.vtk")
    if p["save_stl"]:
        mesh.save(p["output_directory"] / f"city_{prefix}.stl")
        ground_mesh.save(p["output_directory"] / f"ground_{prefix}.stl")
        building_mesh.save(p["output_directory"] / f"building_{prefix}.stl")
    if p["save_obj"]:
        mesh.save(p["output_directory"] / f"city_{prefix}.stl")
        ground_mesh.save(p["output_directory"] / f"ground_{prefix}.stl")
        building_mesh.save(p["output_directory"] / f"building_{prefix}.stl")
