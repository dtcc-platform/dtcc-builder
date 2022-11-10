#!/usr/bin/env python3

import sys

sys.path.append("pybuilder")

from pathlib import Path
import os

from pybuilder.CityModel import CityModel, building_bounds
from pybuilder.PointCloud import PointCloud, calc_las_bounds
from pybuilder.ElevationModel import ElevationModel
from pybuilder.Parameters import load_parameters
from pybuilder.Utils import bounds_intersect, bounds_area

from pybuilder import Meshing


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
        footprint_bounds = building_bounds(
            p["DataDirectory"] / p["BuildingsFileName"], p["DomainMargin"]
        )
        las_bounds = calc_las_bounds(p["DataDirectory"])
        domain_bounds = bounds_intersect(footprint_bounds, las_bounds)
        origin = domain_bounds[:2]
        p["X0"] = origin[0]
        p["Y0"] = origin[1]
        p["XMin"] = 0.0
        p["YMin"] = 0.0
        p["XMax"] = domain_bounds[0] - domain_bounds[2]
        p["YMax"] = domain_bounds[1] - domain_bounds[3]
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
    if p["PointCloudDirectory"] == "":
        p["PointCloudDirectory"] = p["DataDirectory"]
    return p


def generate_citymodel(p):
    origin, domain_bounds = calc_project_domain(p)

    if bounds_area(domain_bounds) < 100:
        print("Domain too small to generate a city model")
        sys.exit(1)

    pc = PointCloud(p["PointCloudDirectory"], domain_bounds)
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
    dsm = ElevationModel(pc, p["ElevationModelResolution"])
    dtm.smooth_elevation_model(p["GroundSmoothing"])

    cm.extract_building_points(pc)

    # 6 is building classification
    if 6 not in pc.used_classifications and p["RANSACOutlierRemover"]:
        cm.building_points_RANSAC_outlier_remover()

    if p["StatisticalOutlierRemover"]:
        cm.building_points_statistical_outlier_remover()

    cm.compute_building_heights(dtm)

    return cm, dtm


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

    if p["WriteVTK"]:
        Meshing.write_VTK_mesh2D(mesh_2D, p["OutputDirectory"] / "Step31Mesh.vtu")

    # Step 3.2: Generate 3D mesh (layer 3D mesh)
    mesh_3D = Meshing.generate_mesh3D(mesh_2D, p["DomainHeight"], p["MeshResolution"])
    if p["WriteVTK"]:
        mesh_3D_boundary = Meshing.extract_boundary3D(mesh_3D)
        Meshing.write_VTK_surface3D(
            mesh_3D_boundary, p["OutputDirectory"] / "Step32Boundary.vtu"
        )
        Meshing.write_VTK_mesh3D(mesh_3D, p["OutputDirectory"] / "Step32Mesh.vtu")

    # Step 3.3: Smooth 3D mesh (set ground height)
    top_height = dtm.mean() + p["DomainHeight"]
    mesh_3D = Meshing.smooth_mesh3D(mesh_3D, cm, dtm, top_height, False)

    if p["WriteVTK"]:
        mesh_3D_boundary = Meshing.extract_boundary3D(mesh_3D)
        Meshing.write_VTK_surface3D(
            mesh_3D_boundary, p["OutputDirectory"] / "Step33Boundary.vtu"
        )
        Meshing.write_VTK_mesh3D(mesh_3D, p["OutputDirectory"] / "Step33Mesh.vtu")

    # Step 3.4: Trim 3D mesh (remove building interiors)
    mesh_3D = Meshing.trim_mesh3D(mesh_3D, mesh_2D, cm, mesh_3D.numLayers)

    # Step 3.5: Smooth 3D mesh (set ground and building heights)
    mesh_3D = Meshing.smooth_mesh3D(mesh_3D, cm, dtm, top_height, True)

    if p["WriteVTK"]:
        mesh_3D_boundary = Meshing.extract_boundary3D(mesh_3D)
        Meshing.write_VTK_surface3D(
            mesh_3D_boundary, p["OutputDirectory"] / "Step35Boundary.vtu"
        )
        Meshing.write_VTK_mesh3D(mesh_3D, p["OutputDirectory"] / "Step35Mesh.vtu")

    return mesh_3D


def main(p, project_path, citymodel_only=False, mesh_only=False):
    # parameters_file, project_path = get_project_paths(args)

    # if not parameters_file.is_file():
    #     print(f"Warning!: cannot find {parameters_file} using default parameters")
    #     p = load_parameters()
    # else:
    #     p = load_parameters(parameters_file)

    p = set_directories(p, project_path)

    if not mesh_only:
        cm, dtm = generate_citymodel(p)
        if p["WriteJSON"]:
            cm.to_JSON(p["OutputDirectory"] / "CityModel.json")
            dtm.to_JSON(p["OutputDirectory"] / "DTM.json", cm.origin)

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

        if p["WriteVTK"]:
            Meshing.write_VTK_mesh3D(volume_mesh, p["OutputDirectory"] / "CityMesh.vtu")
            Meshing.write_VTK_surface3D(
                surface, p["OutputDirectory"] / "CitySurface.vtu"
            )

        if p["WriteSTL"]:
            Meshing.write_surface(surface, p["OutputDirectory"] / "CitySurface.stl")

        if p["WriteOBJ"]:
            Meshing.write_surface(surface, p["OutputDirectory"] / "CitySurface.obj")


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

    args = parser.parse_args()

    parameters_file, project_path = get_project_paths(args.projectpath)

    if not parameters_file.is_file():
        print(f"Warning!: cannot find {parameters_file} using default parameters")
        parameters = load_parameters()
    else:
        parameters = load_parameters(parameters_file)

    main(parameters, project_path, args.citymodel_only, args.mesh_only)
