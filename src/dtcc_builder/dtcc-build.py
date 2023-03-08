#!/usr/bin/env python3

# Copyright (C) 2022 Dag Wästberg
# Licensed under the MIT License
#
# Modified by Anders Logg 2022

import sys

sys.path.append("pybuilder")

from pathlib import Path
import os
import re

import dtcc_io as io
import dtcc_builder as builder
from dtcc_builder import CityModel
from dtcc_builder import ElevationModel
from dtcc_builder import PointCloud
from dtcc_builder import _pybuilder
from dtcc_builder import load_parameters
from dtcc_builder import Meshing

def camel_to_kebab(name):
    name = re.sub("(.)([A-Z][a-z]+)", r"\1-\2", name)
    return re.sub("([a-z0-9])([A-Z])", r"\1-\2", name).lower()


def create_parameters_options(parser):
    type_map = {"str": str, "float": float, "int": int, "bool": bool}

    name_translator = {}
    p = load_parameters()
    for k, v in p.items():
        kebab_name = camel_to_kebab(k)
        name_translator[k] = kebab_name
        name_translator[kebab_name] = k
        val_type = type(v).__name__
        if val_type == "bool":
            parser.add_argument(f"--{kebab_name}", default=v, action="store_true")
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


def calc_project_domain(p):
    if p["AutoDomain"]:
        footprint_bounds = io.citymodel.building_bounds(
            p["DataDirectory"] / p["BuildingsFileName"], p["DomainMargin"]
        )
        las_bounds = io.pointcloud.calc_las_bounds(p["PointCloudDirectory"])
        domain_bounds = io.bounds.bounds_intersect(footprint_bounds, las_bounds)
        origin = domain_bounds[:2]
        p["X0"] = origin[0]
        p["Y0"] = origin[1]
        p["XMin"] = 0.0
        p["YMin"] = 0.0
        p["XMax"] = domain_bounds[2] - domain_bounds[0]
        p["YMax"] = domain_bounds[3] - domain_bounds[1]
    else:
        origin = (p["X0"], p["Y0"])
        domain_bounds = (
            p["X0"] + p["XMin"],
            p["Y0"] + p["YMin"],
            p["X0"] + p["XMax"],
            p["Y0"] + p["YMax"],
        )
    return (origin, domain_bounds)


def set_directories(p, project_path):
    if p["DataDirectory"] == "":
        p["DataDirectory"] = project_path
    if p["OutputDirectory"] == "":
        p["OutputDirectory"] = project_path

    p["DataDirectory"] = Path(p["DataDirectory"])
    p["OutputDirectory"] = Path(p["OutputDirectory"])
    p["OutputDirectory"].mkdir(parents=True,exist_ok=True)
    if p["PointCloudDirectory"] == "":
        p["PointCloudDirectory"] = p["DataDirectory"]
    else:
        p["PointCloudDirectory"] = Path(p["PointCloudDirectory"])
        if not p["PointCloudDirectory"].is_absolute():
            p["PointCloudDirectory"] = p["DataDirectory"] / p["PointCloudDirectory"]

    return p


def generate_citymodel(p):
    origin, domain_bounds = calc_project_domain(p)

    if io.bounds.bounds_area(domain_bounds) < 100:
        print("Domain too small to generate a city model")
        sys.exit(1)

    pc = builder.PointCloud(pointcloud_path = p["PointCloudDirectory"], bounds= domain_bounds)
    if len(pc) == 0:
        print("Error: Point cloud in domain is empty")
        sys.exit(1)
    pc.set_origin(origin)

    pc.remove_global_outliers(p["OutlierMargin"])

    if p["NaiveVegitationFilter"]:
        pc.vegetation_filter()

    cm = CityModel(p["DataDirectory"] / p["BuildingsFileName"])
    cm.set_origin(origin)
    cm.clean_citymodel()

    dtm = ElevationModel(pc, p["ElevationModelResolution"], [2, 9])
    #dsm = ElevationModel(pc, p["ElevationModelResolution"])
    dtm.smooth_elevation_model(p["GroundSmoothing"])

    cm.extract_building_points(pc)

    # 6 is building classification
    # if buildings are classified then we generally don't need this step
    if 6 not in pc.used_classifications and p["RANSACOutlierRemover"]:
        cm.building_points_RANSAC_outlier_remover()

    if p["StatisticalOutlierRemover"]:
        cm.building_points_statistical_outlier_remover()

    cm.compute_building_heights(dtm)

    return cm, pc, dtm


def generate_surface_mesh(cm: CityModel, dtm: ElevationModel, p: dict):
    surfaces = Meshing.generate_surface3D(cm, dtm, p["MeshResolution"])
    ground_surface = surfaces[0]
    building_surfaces = surfaces[1:]
    building_surfaces = Meshing.merge_surfaces3D(building_surfaces)

    return ground_surface, building_surfaces


def generate_volume_mesh(cm: CityModel, dtm: ElevationModel, p: dict):
    if not cm.cleaned:
        cm.clean_citymodel()
    if not cm.simplified:
        print(dtm.bounds)
        cm.simplify_citymodel(dtm.bounds)
    if not cm.calculated_heights:
        cm.compute_building_heights(dtm)

    # Step 3.1: Generate 2D mesh
    mesh_2D = Meshing.generate_mesh2D(cm, dtm.bounds, p["MeshResolution"])

    if p["Debug"] and p["WriteVTK"]:
        Meshing.write_VTK_mesh2D(mesh_2D, p["OutputDirectory"] / "Step31Mesh.vtu")

    # Step 3.2: Generate 3D mesh (layer 3D mesh)
    mesh_3D = Meshing.generate_mesh3D(mesh_2D, p["DomainHeight"], p["MeshResolution"])
    if p["Debug"] and p["WriteVTK"]:
        mesh_3D_boundary = Meshing.extract_boundary3D(mesh_3D)
        Meshing.write_surface(
            mesh_3D_boundary, p["OutputDirectory"] / "Step32Boundary.vtu"
        )
        Meshing.write_volume_mesh(mesh_3D, p["OutputDirectory"] / "Step32Mesh.vtu")

    # Step 3.3: Smooth 3D mesh (set ground height)
    top_height = dtm.mean() + p["DomainHeight"]
    mesh_3D = Meshing.smooth_mesh3D(mesh_3D, cm, dtm, top_height, False)

    if p["Debug"] and p["WriteVTK"]:
        mesh_3D_boundary = Meshing.extract_boundary3D(mesh_3D)
        Meshing.write_VTK_surface3D(
            mesh_3D_boundary, p["OutputDirectory"] / "Step33Boundary.vtu"
        )
        Meshing.write_volume_mesh(mesh_3D, p["OutputDirectory"] / "Step33Mesh.vtu")

    # Step 3.4: Trim 3D mesh (remove building interiors)
    mesh_3D = Meshing.trim_mesh3D(mesh_3D, mesh_2D, cm, mesh_3D.numLayers)

    # Step 3.5: Smooth 3D mesh (set ground and building heights)
    mesh_3D = Meshing.smooth_mesh3D(mesh_3D, cm, dtm, top_height, True)

    if p["Debug"] and p["WriteVTK"]:
        mesh_3D_boundary = Meshing.extract_boundary3D(mesh_3D)
        Meshing.write_VTK_surface3D(
            mesh_3D_boundary, p["OutputDirectory"] / "Step35Boundary.vtu"
        )
        Meshing.write_volume_mesh(mesh_3D, p["OutputDirectory"] / "Step35Mesh.vtu")

    return mesh_3D


def run(p, project_path, citymodel_only=False, mesh_only=False):
    # parameters_file, project_path = get_project_paths(args)

    # if not parameters_file.is_file():
    #     print(f'Warning!: cannot find {parameters_file} using default parameters')
    #     p = load_parameters()
    # else:
    #     p = load_parameters(parameters_file)

    p = set_directories(p, project_path)
    print(p)
    if not mesh_only:

        # Generate city model
        cm, pc, dtm = generate_citymodel(p)

        # Write data to files
        if p["WriteJSON"]:
            cm.to_JSON(p["OutputDirectory"] / "CityModel.json")
            #dtm.to_JSON(p["OutputDirectory"] / "DTM.json", cm.origin)
        if p["WriteProtobuf"]:
            with open(p["OutputDirectory"] / "CityModel.pb", "wb") as dst:
                dst.write(cm.to_protobuf())
            with open(p["OutputDirectory"] / "PointCloud.pb", "wb") as dst:
                dst.write(pc.to_protobuf())

    if not citymodel_only:

        if mesh_only:
            cm = CityModel()
            cm.from_JSON((p["OutputDirectory"] / "CityModel.json"))
            dtm = ElevationModel()
            dtm.from_JSON(p["OutputDirectory"] / "DTM.json")

        ground_surface, builing_surface = generate_surface_mesh(cm, dtm, p)
        volume_mesh = generate_volume_mesh(cm, dtm, p)

        boundary = Meshing.extract_boundary3D(volume_mesh)
        surface = Meshing.extract_open_surface3D(boundary)

        # Write data to files
        if p["WriteVTK"]:
            Meshing.write_volume_mesh(volume_mesh, p["OutputDirectory"] / "CityMesh.vtu")
            Meshing.write_surface(
                surface, p["OutputDirectory"] / "CitySurface.vtu"
            )
        if p["WriteSTL"]:
            Meshing.write_surface(surface, p["OutputDirectory"] / "CitySurface.stl")
        if p["WriteOBJ"]:
            Meshing.write_surface(surface, p["OutputDirectory"] / "CitySurface.obj")
        if p["WriteProtobuf"]:
            Meshing.write_Protobuf_surface3D(
                surface, p["OutputDirectory"] / "CitySurface.pb"
            )
            Meshing.write_Protobuf_surface3D(
                ground_surface, p["OutputDirectory"] / "GroundSurface.pb"
            )


if __name__ == "__main__":
    import sys
    import argparse

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
        print(f"Warning!: cannot find {parameters_file} using default parameters")
        parameters = builder.Parameters.load_parameters()
    else:
        parameters = builder.Parameters.load_parameters(parameters_file)

    parameters = update_parameters_from_options(parameters, args, name_translator)
    run(parameters, project_path, args.citymodel_only, args.mesh_only)
