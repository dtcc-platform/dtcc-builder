# Copyright (C) 2022 Dag Wästberg
# Licensed under the MIT License
#
# Modified by Anders Logg 2022

import sys
import argparse


from pathlib import Path
import os
import re

from google.protobuf.json_format import MessageToJson

import dtcc_io as io
import dtcc_builder as builder

def camel_to_kebab(name):
    name = re.sub("(.)([A-Z][a-z]+)", r"\1-\2", name)
    return re.sub("([a-z0-9])([A-Z])", r"\1-\2", name).lower()


def create_parameters_options(parser):
    type_map = {"str": str, "float": float, "int": int, "bool": bool}

    name_translator = {}
    p = builder.load_parameters()
    for k, v in p.items():
        kebab_name = camel_to_kebab(k)
        name_translator[k] = kebab_name
        name_translator[kebab_name] = k
        val_type = type(v).__name__
        if val_type == "bool":
            parser.add_argument(f"--{kebab_name}",
                                default=v, action="store_true")
        else:
            parser.add_argument(
                f"--{kebab_name}", default=None, type=type_map[val_type]
            )
    return parser, name_translator


def update_parameters_from_options(p, args, name_translator):
    for k, v in p.items():
        print(k)
        snake_name = name_translator.get(k)
        if snake_name is None:
            continue
        snake_name = name_translator[k].replace("-", "_")
        parser_val = getattr(args, snake_name)
        if parser_val is not None:
            p[k] = parser_val
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

    print(p)
    #cm, dtm = None, None
    cm, dtm = builder.create_citymodel(p["DataDirectory"] / p["BuildingsFileName"],p["PointCloudDirectory"],p )
    if p["WriteJSON"]:
        with open(p["OutputDirectory"]/ "CityModel.json", "w") as dst:
            dst.write(MessageToJson(cm))
    if p["WriteProtobuf"]:
        with open(p["OutputDirectory"]/ "CityModel.pb", "wb") as dst:
            dst.write(cm.SerializeToString())
    io.citymodel.write(cm, p["OutputDirectory"] / "CityModel.shp",)
    if not citymodel_only:
        volume_mesh, surface_mesh = builder.create_mesh(p["DataDirectory"] / p["BuildingsFileName"],p["PointCloudDirectory"],p )

        if p["WriteJSON"]:
            with open(p["OutputDirectory"]/ "CitySurface.json", "w"):
                MessageToJson(surface_mesh)
            with open(p["OutputDirectory"]/ "CityMesh.json", "w"):
                MessageToJson(volume_mesh)
        
        if p["WriteProtobuf"]:
            with open(p["OutputDirectory"]/ "CitySurface.pb", "wb") as dst:
                dst.write(surface_mesh.SerilaizeToString())
            with open(p["OutputDirectory"]/ "CityMesh.pb", "wb") as dst:
                dst.write(volume_mesh.SerializeToString())
        if p["WriteVTK"]:
            io.mesh.write(surface_mesh,p["OutputDirectory"]/ "CitySurface.vtk")
            io.mesh.write(volume_mesh,p["OutputDirectory"]/ "CityMesh.vtk", volume_mesh=True)
        if p["WriteOBJ"]:
            io.mesh.write(surface_mesh,p["OutputDirectory"]/ "CitySurface.obj")
            io.mesh.write(volume_mesh,p["OutputDirectory"]/ "CityMesh.obj", volume_mesh=True)
        if p["WriteSTL"]:
            io.mesh.write(surface_mesh,p["OutputDirectory"]/ "CitySurface.stl")
            io.mesh.write(volume_mesh,p["OutputDirectory"]/ "CityMesh.stl", volume_mesh=True)
    
def main():
    parser = argparse.ArgumentParser(
        prog="pybuilder",
        description="Build LoD1 CItyModel mesh fromm footprint and pointcloud",
    )

    parser.add_argument("--citymodel-only", action="store_true")
    parser.add_argument("--mesh-only", action="store_true")
    
    parser.add_argument("projectpath", nargs="?", default=os.getcwd())

    parser, name_translator = create_parameters_options(parser)

    args = parser.parse_args()

    parameters_file, project_path = get_project_paths(args.projectpath)

    if not parameters_file.is_file():
        print(
            f"Warning!: cannot find {parameters_file} using default parameters")
        parameters = builder.Parameters.load_parameters(None, project_path)
    else:
        parameters = builder.Parameters.load_parameters(parameters_file, project_path)

    parameters = update_parameters_from_options(
        parameters, args, name_translator)
    run(parameters, args.citymodel_only, args.mesh_only)
