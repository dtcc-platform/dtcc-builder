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


def set_directories(p, project_path):
    if p["data_directory"] == "":
        p["data_directory"] = project_path
    if p["output_directory"] == "":
        p["output_directory"] = p["data_directory"]

    p["data_directory"] = Path(p["data_directory"])
    p["output_directory"] = Path(p["output_directory"])
    p["output_directory"].mkdir(parents=True, exist_ok=True)
    if p["pointcloud_directory"] == "":
        p["pointcloud_directory"] = p["data_directory"]
    else:
        p["pointcloud_directory"] = Path(p["pointcloud_directory"])
        if not p["pointcloud_directory"].is_absolute():
            p["pointcloud_directory"] = p["data_directory"] / p["pointcloud_directory"]

    return p


def load_parameters(file_path=None, project_path="."):
    if file_path is None:
        p = builder.parameters.default()

    else:
        file_path = Path(file_path)
        print(file_path)
        if file_path.is_dir():
            file_path = file_path / "Parameters.json"
        if not file_path.exists():
            print(
                f"Parameters file {file_path} does not exist, using default parameters"
            )
            p = builder.parameters.default()
        else:
            with open(file_path) as src:
                loaded_parameters = json.load(src)
            p = builder.parameters.default()
            p.update(loaded_parameters)
    p = set_directories(p, project_path)
    return p


def create_parameters_options(parser):
    type_map = {
        "str": str,
        "float": float,
        "int": int,
        "bool": bool,
    }

    p = load_parameters()
    for key, value in p.items():
        val_type = type(value).__name__
        if isinstance(value, bool):
            parser.add_argument(
                f"--{key}",
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
            parser.add_argument(
                f"--{key}", default=None, type=type_map.get(val_type, str)
            )
    return parser


def update_parameters_from_options(parameters, args):
    for key, value in parameters.items():
        parser_value = getattr(args, key)
        if parser_value is not None:
            parameters[key] = parser_value
    return parameters


def print_parameters(p):
    "Pretty-print parameters"
    info("Printing parameters")
    keys = sorted(p.keys())
    n = max([len(key) for key in keys])
    for key in sorted(p.keys()):
        print(f"  {key}: {' '*(n - len(key) - 1)} {p[key]}")


def get_project_paths(path):
    if not path:
        project_path = Path.cwd()
        parameters_file = project_path / PARAMETERS_FILE
    else:
        arg_path = Path(path).resolve()
        if not arg_path.exists():
            raise FileNotFoundError(f"cannot find project path {arg_path}")
        if arg_path.is_file():
            project_path = arg_path.parent
            parameters_file = arg_path
        else:
            project_path = arg_path
            parameters_file = project_path / PARAMETERS_FILE

    return (parameters_file, project_path)


def run(p, citymodel_only, mesh_only):
    building_file = p["data_directory"] / p["buildings_filename"]
    if not building_file.exists():
        raise FileNotFoundError(f"cannot find building file {building_file}")
    pointcloud_file = p["pointcloud_directory"]
    if not pointcloud_file.exists():
        raise FileNotFoundError(f"cannot find point cloud file {pointcloud_file}")
    info(
        f"creating citymodel from Building file: {building_file} and pointcloud {pointcloud_file}"
    )

    origin, project_bounds = builder.build.calculate_project_domain(
        building_file, pointcloud_file, p
    )

    info(f"project bounds: {project_bounds.tuple}")

    cm = io.load_citymodel(
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

    dem_raster = builder.build.generate_dem(
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

    io.save_citymodel(
        cm,
        p["output_directory"] / "CityModel.shp",
    )

    if p["write_protobuf"]:
        io.save_citymodel(
            cm,
            p["output_directory"] / "CityModel.pb",
        )

    if p["write_json"]:
        io.save_citymodel(
            cm,
            p["output_directory"] / "CityModel.json",
        )

    if not citymodel_only:
        pass

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
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        prog="dtcc-builder",
        description="Build LOD1 city mesh(es) from building footprints and pointcloud",
    )

    parser.add_argument("--citymodel-only", action="store_true")
    parser.add_argument("--mesh-only", action="store_true")
    parser.add_argument("projectpath", nargs="?", default=os.getcwd())

    parser = create_parameters_options(parser)

    args = parser.parse_args()

    parameters_file, project_path = get_project_paths(args.projectpath)

    if not parameters_file.is_file():
        warning(f"Unable to load {parameters_file}; using default parameters")
        parameters = load_parameters(None, project_path)
    else:
        info(f"Loading parameters from {parameters_file}")
        parameters = load_parameters(parameters_file, project_path)

    parameters = update_parameters_from_options(parameters, args)

    # Pretty-print parameters
    print_parameters(parameters)

    # Run the builder
    run(parameters, args.citymodel_only, args.mesh_only)
