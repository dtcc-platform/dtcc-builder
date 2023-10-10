# Copyright (C) 2022 Dag Wästberg
# Licensed under the MIT License
#
# Modified by Anders Logg 2023
#
# This module provides the main methods of DTCC Builder.

import numpy as np
import dtcc_wrangler
from pathlib import Path
from typing import Tuple, List
from time import time


import dtcc_model as model
import dtcc_io as io
from .logging import info, warning, error
from . import _dtcc_builder
from . import city_methods
from . import model as builder_model
from . import parameters as builder_parameters
from . import meshing


def calculate_bounds(
    buildings_path, pointcloud_path, parameters: dict = None
) -> Tuple[Tuple[float, float], model.Bounds]:
    """
    Compute the domain bounds based on building footprints and pointcloud data.

    Parameters
    ----------
    `buildings_path` : str
        Path to the building footprints file.
    `pointcloud_path` : str
        Path to the point cloud data.
    `parameters` : dict, optional
        A dictionary of parameters for the computation, by default None.

    Returns
    -------
    `origin` : Tuple[float, float]
        Tuple containing the origin coordinates.
    `bounds` : dtcc_model.Bounds
    Tuple[Tuple[float, float], dtcc_model.Bounds]
        The computed domain bounds.

    """

    # Get parameters
    p = parameters or builder_parameters.default()

    # Compute domain bounds automatically or from parameters
    if p["auto_domain"]:
        info("Calculating domain bounds automatically...")
        city_bounds = io.city.building_bounds(buildings_path, p["domain_margin"])
        info(f"Footprint bounds: {city_bounds}")
        pointcloud_bounds = io.pointcloud.calc_las_bounds(pointcloud_path)
        info(f"Point cloud bounds: {pointcloud_bounds}")
        city_bounds.intersect(pointcloud_bounds)
        info(f"Intersected bounds: {city_bounds}")
        bounds = city_bounds.tuple
        origin = bounds[:2]
        p["x0"] = origin[0]
        p["y0"] = origin[1]
        p["x_min"] = 0.0
        p["y_min"] = 0.0
        p["x_max"] = bounds[2] - bounds[0]
        p["y_max"] = bounds[3] - bounds[1]
    else:
        info("Calculating domain bounds from parameters...")
        origin = (p["x0"], p["y0"])
        bounds = (
            p["x0"] + p["x_min"],
            p["y0"] + p["y_min"],
            p["x0"] + p["x_max"],
            p["y0"] + p["y_max"],
        )

    # Set bounds
    bounds = model.Bounds(
        xmin=bounds[0],
        ymin=bounds[1],
        xmax=bounds[2],
        ymax=bounds[3],
    )

    return origin, bounds


def build_city(
    city: model.City,
    point_cloud: model.PointCloud,
    bounds: model.Bounds,
    parameters: dict = None,
) -> model.City:
    """
    Build a city model from building footprints and point cloud data.

    Parameters
    ----------
    `city` : dtcc_model.City
        The city model to be built.
    `point_cloud` : dtcc_model.PointCloud
        The point cloud data associated with the city.
    `bounds` : dtcc_model.Bounds
        The bounds of the city domain.
    `parameters` : dict, optional
        A dictionary of parameters for the computation, by default None.

    Returns
    -------
    dtcc_model.City
        The constructed city model.

    Developer notes
    --------------
    Consider introducing a new class named Footprints
    so that a city can be built from footprints and point cloud data.
    It is somewhat strange that the input to this function is a city.
    """

    info("Building city...")

    # Get parameters
    p = parameters or builder_parameters.default()

    start_time = time()
    # Remove outliers from point cloud
    point_cloud = point_cloud.remove_global_outliers(p["outlier_margin"])
    # FIXME: Don't modify incoming data (city)
    print(f"BBB: Removing outliers took {time() - start_time} seconds")

    # FIXME: Why are we not calling clean_city?
    # Should be callled with min_vertex_distance/2.

    start_time = time()
    # Build elevation model
    city = city.terrain_from_pointcloud(
        point_cloud,
        p["elevation_model_resolution"],
        p["elevation_model_window_size"],
        ground_only=True,
    )
    print(f"BBB: Building elevation model took {time() - start_time} seconds")

    start_time = time()
    city = city.simplify_buildings(p["min_building_distance"] / 2)
    print(f"BBB: Simplifying buildings took {time() - start_time} seconds")

    # Compute building points
    start_time = time()
    city = city_methods.compute_building_points(
        city,
        point_cloud,
        p["ground_margin"],
        p["outlier_margin"],
        p["statistical_outlier_remover"],
        p["roof_outlier_neighbors"],
        p["roof_outlier_margin"],
        p["ransac_outlier_remover"],
        p["ransac_outlier_margin"],
        p["ransac_iterations"],
    )
    print(f"BBB: Computing building points took {time() - start_time} seconds")

    # Compute building heights
    start_time = time()
    city = city_methods.compute_building_heights(
        city, p["min_building_height"], p["roof_percentile"]
    )
    print(f"BBB: Computing building heights took {time() - start_time} seconds")

    return city


def build_terrain_mesh(city: model.City, parameters: dict = None):
    """
    Build terrain mesh for city.

    This function builds a terrain mesh for the city

    Parameters
    ----------
    `city` : dtcc_model.City
        The city model for which to generate meshes.
    `parameters` : dict, optional
        A dictionary of parameters for the computation, by default None.

    Returns
    -------
    `dtcc_ground_mesh` : dtcc_model.Mesh
        The ground mesh
    """
    info("Building terrain mesh...")

    # Get parameters
    p = parameters or builder_parameters.default()

    simple_city = city.remove_small_buildings(
        p["min_building_size"]
    ).simplify_buildings(p["min_building_distance"] / 2)

    return meshing.terrain_mesh(
        simple_city, p["mesh_resolution"], p["ground_smoothing"]
    )


def build_building_meshes(
    city: model.City, parameters: dict = None
) -> List[model.Mesh]:
    """
    Build meshes for each building in the city.

    Returns a list of Meshes, one for each building in the city.

    Parameters
    ----------
    `city` : dtcc_model.City
        The city model for which to generate meshes.
    `parameters` : dict, optional
        A dictionary of parameters for the computation, by default None.

    Returns
    -------

    `dtcc_building_meshes` : [dtcc_model.Mesh]
        Meshes of all the buildings

    """
    info("meshing buildings...")

    # Get parameters
    p = parameters or builder_parameters.default()

    simple_city = city.remove_small_buildings(
        p["min_building_size"]
    ).simplify_buildings(p["min_building_distance"] / 2)

    building_meshes = meshing.building_meshes(
        city, p["mesh_resolution"], cap_base=True, per_floor=True
    )

    return building_meshes


def build_volume_mesh(
    city: model.City, parameters: dict = None
) -> Tuple[model.VolumeMesh, model.Mesh]:
    """
    Build volume mesh for city.

    This function builds a boundary-conforming volume mesh for the city,
    returning both the volume mesh and its corresponding boundary mesh.

    Parameters
    ----------
    `city` : dtcc_model.City
        The city model for which to generate the volume mesh.
    `parameters` : dict, optional
        A dictionary of parameters for the computation, by default None.

    Returns
    -------
    `volume_mesh` : dtcc_model.VolumeMesh
        The city's volume mesh
    `volume_mesh_boundary` : dtcc_model.Mesh
        The city's boundary mesh (boundary of the volume mesh).

    """
    info("Building volume mesh for city...")

    # Get parameters
    p = parameters or builder_parameters.default()
    simple_city = city.merge_buildings(
        p["min_building_distance"]
    ).remove_small_buildings(p["min_building_size"])

    # Convert to builder model
    simple_builder_city = builder_model.create_builder_city(simple_city)
    builder_dem = builder_model.raster_to_builder_gridfield(city.terrain)

    # Step 3.1: Build ground mesh
    ground_mesh = _dtcc_builder.build_ground_mesh(
        simple_builder_city,
        city.bounds.xmin,
        city.bounds.ymin,
        city.bounds.xmax,
        city.bounds.ymax,
        p["mesh_resolution"],
        True,
    )
    _debug(ground_mesh, "3.1", p)

    # Step 3.2: Layer ground mesh
    volume_mesh = _dtcc_builder.layer_ground_mesh(
        ground_mesh, p["domain_height"], p["mesh_resolution"]
    )
    _debug(volume_mesh, "3.2", p)

    # Step 3.3: Smooth volume mesh (set ground height)
    top_height = p["domain_height"] + city.terrain.data.mean()
    volume_mesh = _dtcc_builder.smooth_volume_mesh(
        volume_mesh,
        simple_builder_city,
        builder_dem,
        top_height,
        False,
        p["smoothing_max_iterations"],
        p["smoothing_relative_tolerance"],
    )
    _debug(volume_mesh, "3.3", p)

    # Step 3.4: Trim volume mesh (remove building interiors)
    volume_mesh = _dtcc_builder.trim_volume_mesh(
        volume_mesh, ground_mesh, simple_builder_city
    )
    _debug(volume_mesh, "3.4", p)

    # Step 3.5: Smooth volume mesh (set ground and building heights)
    volume_mesh = _dtcc_builder.smooth_volume_mesh(
        volume_mesh,
        simple_builder_city,
        builder_dem,
        top_height,
        True,
        p["smoothing_max_iterations"],
        p["smoothing_relative_tolerance"],
    )
    _debug(volume_mesh, "3.5", p)

    # Compute boundary mesh
    volume_mesh_boundary = _dtcc_builder.compute_boundary_mesh(volume_mesh)

    # Convert back to DTCC model
    dtcc_volume_mesh = builder_model.builder_volume_mesh_to_volume_mesh(volume_mesh)
    dtcc_volume_mesh_boundary = builder_model.builder_mesh_to_mesh(volume_mesh_boundary)

    return dtcc_volume_mesh, dtcc_volume_mesh_boundary


def build_city_surface_mesh(
    city: model.City, parameters: dict = None, merge_meshes=True
) -> model.Mesh:
    """
    Build surface mesh for city.
    Parameters
    ----------
    `city` : dtcc_model.City
        The city model for which to generate the volume mesh.
    `parameters` : dict, optional
        A dictionary of parameters for the computation, by default None.
    Returns
    -------
    `surface_mesh` : dtcc_model.Mesh
        The city's surface mesh
    """
    p = parameters or builder_parameters.default()
    mesh_resolutiom = p["mesh_resolution"]
    ground_smoothing = p["ground_smoothing"]
    surface_mesh = meshing.city_surface_mesh(
        city, mesh_resolutiom, ground_smoothing, merge_meshes
    )
    return surface_mesh


def build(parameters: dict = None) -> None:
    """
    Build city and city meshes.

    This function reads data from the specified data directory and builds a city and its corresponding meshes.
    The same thing can be accomplished by calling the individual build_*
    functions, but this function is provided as a convenience
    and takes care of loading and saving data to files.

    Parameters
    ----------
    `parameters` : dict, optional
        A dictionary of parameters for the computation, by default None.

    """

    # Get parameters
    p = parameters or builder_parameters.default()

    # Get paths
    if p["data_directory"]:
        data_directory = Path(p["data_directory"])
    else:
        data_directory = Path.cwd()
    if p["output_directory"]:
        output_directory = Path(p["output_directory"])
    else:
        output_directory = data_directory
    if p["pointcloud_directory"]:
        pointcloud_path = Path(p["pointcloud_directory"])
    else:
        pointcloud_path = data_directory
    buildings_path = data_directory / p["buildings_filename"]
    if not buildings_path.exists():
        error(f"Unable to read buildings file {buildings_path}")
    if not pointcloud_path.exists():
        error(f"Unable to read pointcloud directory {pointcloud_path}")

    # Compute domain bounds
    origin, bounds = calculate_bounds(buildings_path, pointcloud_path, p)
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
    city = build_city(city, point_cloud, bounds, p)

    # Save city to file
    city_name = output_directory / "city"
    if p["save_protobuf"]:
        io.save_city(city, city_name.with_suffix(".pb"))
    if p["save_shp"]:
        io.save_city(city, city_name.with_suffix(".shp"))
    if p["save_json"]:
        io.save_city(city, city_name.with_suffix(".json"))

    # Build mesh
    if p["build_mesh"]:
        ground_mesh = build_terrain_mesh(city, p)
        building_meshes = build_building_meshes(city, p)
        building_mesh = meshing.merge_meshes(building_meshes)
        surface_mesh = build_city_surface_mesh(city, p)

        # Save meshes to file
        ground_mesh_name = output_directory / "ground_mesh"
        building_mesh_name = output_directory / "building_mesh"
        surface_mesh_name = output_directory / "surface_mesh"
        if p["save_protobuf"]:
            ground_mesh.save(ground_mesh_name.with_suffix(".pb"))
            building_mesh.save(building_mesh_name.with_suffix(".pb"))
            surface_mesh.save(surface_mesh_name.with_suffix(".pb"))
        if p["save_vtk"]:
            ground_mesh.save(ground_mesh_name.with_suffix(".vtu"))
            building_mesh.save(building_mesh_name.with_suffix(".vtu"))
            surface_mesh.save(surface_mesh_name.with_suffix(".vtu"))
        if p["save_stl"]:
            ground_mesh.save(ground_mesh_name.with_suffix(".stl"))
            building_mesh.save(building_mesh_name.with_suffix(".stl"))
            surface_mesh.save(surface_mesh_name.with_suffix(".stl"))
        if p["save_obj"]:
            ground_mesh.save(ground_mesh_name.with_suffix(".obj"))
            building_mesh.save(building_mesh_name.with_suffix(".obj"))
            surface_mesh.save(surface_mesh_name.with_suffix(".obj"))
    # Build volume mesh
    if p["build_volume_mesh"]:
        volume_mesh, volume_mesh_boundary = build_volume_mesh(city, p)

        # Save meshes to file
        volume_mesh_name = output_directory / "volume_mesh"
        volume_mesh_boundary_name = output_directory / "volume_mesh_boundary"
        if p["save_protobuf"]:
            volume_mesh.save(volume_mesh_name.with_suffix(".pb"))
            volume_mesh_boundary.save(volume_mesh_boundary_name.with_suffix(".pb"))
        if p["save_vtk"]:
            volume_mesh.save(volume_mesh_name.with_suffix(".vtu"))
            volume_mesh_boundary.save(volume_mesh_boundary_name.with_suffix(".vtu"))
        if p["save_stl"]:
            volume_mesh_boundary.save(volume_mesh_boundary_name.with_suffix(".stl"))
        if p["save_obj"]:
            volume_mesh_boundary.save(volume_mesh_boundary_name.with_suffix(".obj"))


def _debug(mesh, step, p):
    "Debug volume meshing"
    if not p["debug"]:
        return
    if isinstance(mesh, _dtcc_builder.Mesh):
        mesh = builder_model.builder_mesh_to_mesh(mesh)
    else:
        mesh = builder_model.builder_volume_mesh_to_volume_mesh(mesh)
    mesh.save(p["output_directory"] / f"mesh_step{step}.vtu")
