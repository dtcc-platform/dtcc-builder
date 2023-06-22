# Copyright (C) 2022 Dag Wästberg
# Licensed under the MIT License
#
# Modified by Anders Logg 2022
#
# This is the main command line script provided by DTCC Builder.
# It builds city models and/or meshes from raw data in a a given
# project directory.


import sys, os, re, json, argparse
from pathlib import Path


import dtcc_io as io
import dtcc_builder as builder
from dtcc_common import info, warning

PARAMETERS_FILE = "parameters.json"


def parse_command_line():
    "Parse command line arguments"

    # Create parser
    parser = argparse.ArgumentParser(
        prog="dtcc-builder",
        description="Build LOD1 city mesh(es) from building footprints and pointcloud",
    )

    # Add main arguments
    parser.add_argument("--city-only", action="store_true")
    parser.add_argument("--mesh-only", action="store_true")
    parser.add_argument("path", nargs="?", default=os.getcwd())

    # Add parameter arguments (note special handling of booleans)
    for key, value in builder.parameters.default().items():
        if isinstance(value, bool):
            parser.add_argument(
                f"--{key}",
                dest=key,
                action="store_true",
                default=None,
                help=f"Turn on {key}",
            )
            parser.add_argument(
                f"--no-{key}",
                dest=key,
                default=None,
                action="store_false",
                help=f"Turn off {key}",
            )
        else:
            parser.add_argument(f"--{key}", default=None, type=type(value))

    return parser.parse_args()


def load_parameters(args):
    "Load parameters and override by command line arguments"

    # Load parameters from or directory if provided, else use default parameters
    parameters = builder.parameters.default()
    if args.path is None:
        info("Using default parameters")
    else:
        path = Path(args.path)
        if path.is_dir():
            path = path / PARAMETERS_FILE
        if path.exists():
            info(f"Loading parameters from {path}")
            with open(path) as f:
                _parameters = json.load(f)
                parameters.update(_parameters)
        else:
            warning(f"Unable to load {path}; using default parameters")

    # Override parameters with command line arguments
    for key, value in parameters.items():
        parser_value = getattr(args, key)
        if parser_value is not None:
            info(f"Overriding parameter {key} with argument {parser_value}")
            parameters[key] = parser_value

    # Set parameters for directories
    set_directory_parameters(parameters, args.path)

    # Pretty-print parameters
    info("Printing parameters")
    keys = sorted(parameters.keys())
    n = max([len(key) for key in keys]) if keys else 0
    for key in sorted(keys):
        print(f"  {key}: {' '*(n - len(key) - 1)} {parameters[key]}")

    return parameters


def set_directory_parameters(parameters, path):
    "Set parameters for directories"

    # Shortcut
    p = parameters

    # Set data_directory
    data_directory = p["data_directory"]
    if data_directory == "":
        data_directory = path
    data_directory = Path(data_directory)

    # Set pointcloud_directory
    pointcloud_directory = p["pointcloud_directory"]
    if pointcloud_directory == "":
        pointcloud_directory = data_directory
    pointcloud_directory = Path(pointcloud_directory)

    # Set output_directory
    output_directory = p["output_directory"]
    if output_directory == "":
        output_directory = data_directory
    output_directory = Path(output_directory)

    # Set parameters
    p["data_directory"] = data_directory
    p["pointcloud_directory"] = pointcloud_directory
    p["output_directory"] = output_directory

    return p


def run(p, city_only, mesh_only):
    building_file = p["data_directory"] / p["buildings_filename"]
    if not building_file.exists():
        raise FileNotFoundError(f"cannot find building file {building_file}")
    pointcloud_file = p["pointcloud_directory"]
    if not pointcloud_file.exists():
        raise FileNotFoundError(f"cannot find point cloud file {pointcloud_file}")
    info(
        f"creating city from Building file: {building_file} and pointcloud {pointcloud_file}"
    )

    origin, project_bounds = builder.build.calculate_project_domain(
        building_file, pointcloud_file, p
    )

    info(f"project bounds: {project_bounds.tuple}")

    cm = io.load_city(
        building_file,
        uuid_field=p["uuid_field"],
        height_field=p["height_field"],
        bounds=project_bounds,
    )
    pc = io.load_pointcloud(
        pointcloud_file,
        bounds=project_bounds,
    )
    pc = pc.remove_global_outliers(p["outlier_margin"])

    dem_raster = builder.build.build_dem(
        pc, project_bounds, p["elevation_model_resolution"]
    )
    cm.terrain = dem_raster

    if not ["statistical_outlier_remover"]:
        p["roof_outlier_margin"] = 0

    if not p["ransac_outlier_remover"]:
        p["ransac_iterations"] = 0

    cm = builder.build.extract_buildingpoints(
        cm,
        pc,
        p["ground_margin"],
        p["outlier_margin"],
        p["roof_outlier_margin"],
        p["outlier_neighbors"],
        p["ransac_outlier_margin"],
        p["ransac_iterations"],
    )

    cm = builder.build.calculate_building_heights(
        cm, p["roof_percentile"], p["min_building_height"], overwrite=True
    )

    io.save_city(
        cm,
        p["output_directory"] / "City.shp",
    )

    if p["write_protobuf"]:
        io.save_city(
            cm,
            p["output_directory"] / "City.pb",
        )

    if p["write_json"]:
        io.save_city(
            cm,
            p["output_directory"] / "City.json",
        )

    if not city_only:

        volume_mesh, surface_mesh = builder.build.build_mesh(
            cm,
            p["mesh_resolution"],
            p["domain_height"],
            p["min_building_distance"],
            p["min_vertex_distance"],
            p["debug"],
        )

        ground_surface, buildings = builder.build.build_surface_meshes(
            cm,
            p["min_building_distance"],
            p["min_vertex_distance"],
            p["mesh_resolution"],
        )

        if p["write_protobuf"]:
            surface_mesh.save(p["output_directory"] / "CitySurface.pb")
            ground_surface.save(p["output_directory"] / "GroundSurface.pb")
            buildings.save(p["output_directory"] / "Buildings.pb")

        if p["write_vtk"]:
            surface_mesh.save(p["output_directory"] / "CitySurface.vtk")
            volume_mesh.save(p["output_directory"] / "CityMesh.vtk")

        if p["write_stl"]:
            surface_mesh.save(p["output_directory"] / "CitySurface.stl")
            ground_surface.save(p["output_directory"] / "GroundSurface.stl")
            buildings.save(p["output_directory"] / "Buildings.stl")


def main():
    # Parse command line
    args = parse_command_line()

    # Load parameters
    parameters = load_parameters(args)

    # Run builder
    run(parameters, args.city_only, args.mesh_only)
