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


def get_project_paths(args):
    if len(args) == 1:
        project_path = Path.cwd()
        parameters_file = project_path / "Parameters.json"
    else:
        arg_path = Path(args[1])
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


def main(args):
    parameters_file, project_path = get_project_paths(args)

    if not parameters_file.is_file():
        print(f"Warning!: cannot find {parameters_file} using default parameters")
        p = load_parameters()
    else:
        p = load_parameters(parameters_file)

    if p["DataDirectory"] == "":
        p["DataDirectory"] = project_path
    if p["OutputDirectory"] == "":
        p["OutputDirectory"] = project_path

    p["DataDirectory"] = Path(p["DataDirectory"])
    p["OutputDirectory"] = Path(p["OutputDirectory"])
    if p["PointCloudDirectory"] == "":
        p["PointCloudDirectory"] = p["DataDirectory"]

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
    

    cm = CityModel(p["DataDirectory"] / p["BuildingsFileName"] )
    cm.set_origin(origin)
    cm.clean_citymodel()

    dtm = ElevationModel(pc,p["ElevationModelResolution"],[2,9])
    dsm = ElevationModel(pc,p["ElevationModelResolution"])
    dtm.smooth_elevation_model(p["GroundSmoothing"])

    cm.extract_building_points(pc)

    #6 is building classification
    if 6 not in pc.used_classifications and p["RANSACOutlierRemover"]:
        cm.building_points_RANSAC_outlier_remover()
    
    if p["StatisticalOutlierRemover"]:
        cm.building_points_statistical_outlier_remover()

    cm.compute_building_heights(dtm)

    cm.to_JSON(p["OutputDirectory"] / "CityModel.json")


if __name__ == "__main__":
    import sys

    args = sys.argv
    print(args)
    if "-h" in args:
        print("Usage: dtcc-generate-citymodel.py Parameters.json")
        sys.exit(0)
    main(args)
