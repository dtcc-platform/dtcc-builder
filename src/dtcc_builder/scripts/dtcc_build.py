# Copyright (C) 2022 Dag Wästberg
# Licensed under the MIT License
#
# Modified by Anders Logg 2022

import sys
import argparse


from pathlib import Path
import os
import re
import logging

from google.protobuf.json_format import MessageToJson

import dtcc_io as io
import dtcc_builder as builder


def camel_to_kebab(name):
    name = re.sub("(.)([A-Z][a-z]+)", r"\1-\2", name)
    return re.sub("([a-z0-9])([A-Z])", r"\1-\2", name).lower()


def create_parameters_options(parser):
    type_map = {
        "str": str,
        "float": float,
        "int": int,
        "bool": bool,
    }

    name_translator = {}
    p = builder.load_parameters()
    for k, v in p.items():
        kebab_name = camel_to_kebab(k)
        name_translator[k] = kebab_name
        name_translator[kebab_name] = k
        val_type = type(v).__name__
        if val_type == "bool":
            # print(f"bool arg: {kebab_name}")
            parser.add_argument(
                f"--{kebab_name}",
                action="store_true",
                default=None,
                help=f"Turn on {k}",
            )
            parser.add_argument(
                f"--no-{kebab_name}",
                dest=kebab_name.replace("-", "_"),
                default=None,
                action="store_false",
                help=f"Turn off {k}",
            )
        else:
            parser.add_argument(
                f"--{kebab_name}", default=None, type=type_map.get(val_type, str)
            )
    return parser, name_translator


def update_parameters_from_options(p, args, name_translator):
    for k, v in p.items():
        # print(k)
        snake_name = name_translator.get(k)
        # print(snake_name)
        if snake_name is None:
            continue
        snake_name = name_translator[k].replace("-", "_")
        # print(snake_name)
        parser_val = getattr(args, snake_name)
        # print(parser_val)
        if parser_val is not None:
            p[k] = parser_val
        # print("---")
    return p


def get_project_paths(path):
    if not path:
        project_path = Path.cwd()
        parameters_file = project_path / "Parameters.json"
    else:
        arg_path = Path(path).resolve()
        if not arg_path.exists():
            raise FileNotFoundError(f"cannot find project path {arg_path}")
        if arg_path.is_file():
            project_path = arg_path.parent
            parameters_file = arg_path
        else:
            project_path = arg_path
            parameters_file = project_path / "Parameters.json"

    return (parameters_file, project_path)


def run(p, citymodel_only, mesh_only):
    building_file = p["data-directory"] / p["buildings-filename"]
    if not building_file.exists():
        raise FileNotFoundError(f"cannot find building file {building_file}")
    pointcloud_file = p["pointcloud-directory"]
    if not pointcloud_file.exists():
        raise FileNotFoundError(f"cannot find point cloud file {pointcloud_file}")
    logging.info(
        f"creating citymodel from Building file: {building_file} and pointcloud {pointcloud_file}"
    )

    origin, project_bounds = builder.build.calculate_project_domain(
        building_file, pointcloud_file, p
    )

    print(f"project bounds: {project_bounds.tuple}")

    cm = io.load_citymodel(
        building_file,
        uuid_field=p["uuid-field"],
        height_field=p["height-field"],
        bounds=project_bounds,
    )
    pc = io.load_pointcloud(
        pointcloud_file,
        bounds=project_bounds,
    )

    pc = pc.remove_global_outliers(p["outlier-margin"])

    dem_raster = builder.build.generate_dem(
        pc, project_bounds, p["elevation-model-resolution"]
    )
    cm.terrain = dem_raster

    if not ["statistical-outlier-remover"]:
        p["roof-outlier-margin"] = 0

    if not p["ransac-outlier-remover"]:
        p["ransac-iterations"] = 0

    cm = builder.build.extract_buildingpoints(
        cm,
        pc,
        p["ground-margin"],
        p["outlier-margin"],
        p["roof-outlier-margin"],
        p["outlier-neighbors"],
        p["ransac-outlier-margin"],
        p["ransac-iterations"],
    )

    cm = builder.build.calculate_building_heights(
        cm, p["roof-percentile"], p["min-building-height"], overwrite=True
    )

    io.save_citymodel(
        cm,
        p["output-directory"] / "CityModel.shp",
    )

    if p["write-protobuf"]:
        io.save_citymodel(
            cm,
            p["output-directory"] / "CityModel.pb",
        )

    if p["write-json"]:
        io.save_citymodel(
            cm,
            p["output-directory"] / "CityModel.json",
        )

    if not citymodel_only:
        pass

        volume_mesh, surface_mesh = builder.build.build_mesh(
            cm,
            p["mesh-resolution"],
            p["domain-height"],
            p["min-building-distance"],
            p["min-vertex-distance"],
            p["debug"],
        )

        ground_surface, buildings = builder.build.build_surface_meshes(
            cm,
            p["min-building-distance"],
            p["min-vertex-distance"],
            p["mesh-resolution"],
        )

        if p["write-protobuf"]:
            surface_mesh.save(p["output-directory"] / "CitySurface.pb")
            ground_surface.save(p["output-directory"] / "GroundSurface.pb")
            buildings.save(p["output-directory"] / "Buildings.pb")

        if p["write-vtk"]:
            surface_mesh.save(p["output-directory"] / "CitySurface.vtk")
            volume_mesh.save(p["output-directory"] / "CityMesh.vtk")

        if p["write-stl"]:
            surface_mesh.save(p["output-directory"] / "CitySurface.stl")
            ground_surface.save(p["output-directory"] / "GroundSurface.stl")
            buildings.save(p["output-directory"] / "Buildings.stl")


def main():
    parser = argparse.ArgumentParser(
        prog="dtcc-builder",
        description="Build LoD1 CItyModel mesh fromm footprint and pointcloud",
    )

    parser.add_argument("--citymodel-only", action="store_true")
    parser.add_argument("--mesh-only", action="store_true")

    parser.add_argument("projectpath", nargs="?", default=os.getcwd())

    parser, name_translator = create_parameters_options(parser)

    args = parser.parse_args()

    parameters_file, project_path = get_project_paths(args.projectpath)

    if not parameters_file.is_file():
        print(f"Warning!: cannot find {parameters_file} using default parameters")
        parameters = builder.Parameters.load_parameters(None, project_path)
    else:
        parameters = builder.Parameters.load_parameters(parameters_file, project_path)

    parameters = update_parameters_from_options(parameters, args, name_translator)
    print(parameters)
    run(parameters, args.citymodel_only, args.mesh_only)
