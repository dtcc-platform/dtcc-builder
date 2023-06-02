import dtcc_builder as builder
import dtcc_builder.build
import dtcc_builder.builder_datamodel as builder_datamodel
import dtcc_model as model
from dtcc_model import dtcc_pb2 as proto
from dtcc_builder import _pybuilder
import dtcc_io as io
from pathlib import Path
import numpy as np
from affine import Affine
from typing import Union, Tuple

import rasterio.transform

from pypoints2grid import points2grid

from time import time


def calculate_project_domain(params_or_footprint, las_path=None, params=None):
    if las_path is None and params is None:
        p = params_or_footprint
        building_footprint_path = p["DataDirectory"] / p["BuildingsFileName"]
        pointcloud_path = p["PointCloudDirectory"]
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
    if p["AutoDomain"]:
        footprint_bounds = io.citymodel.building_bounds(
            building_footprint_path, p["DomainMargin"]
        )
        las_bounds = io.pointcloud.calc_las_bounds(pointcloud_path)
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
    domain_bounds = model.Bounds(
        xmin=domain_bounds[0],
        ymin=domain_bounds[1],
        xmax=domain_bounds[2],
        ymax=domain_bounds[3],
    )
    return (origin, domain_bounds)


def generate_dem(
    pointcloud: model.PointCloud, bounds, cell_size: float, window_size: int = 3
) -> model.Raster:
    if (
        len(pointcloud.classification) == len(pointcloud.points)
    ) and 2 in pointcloud.used_classifications():
        ground_point_idx = np.where(np.isin(pointcloud.classification, [2, 9]))[0]
        ground_points = pointcloud.points[ground_point_idx]
    else:
        ground_points = pointcloud.points
    print(f"generating dem with bounds {bounds.tuple}")
    dem = points2grid(ground_points, cell_size, bounds.tuple, window_size=window_size)

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


def extract_buildingpoints(
    citymodel,
    pointcloud,
    ground_padding=2.0,
    ground_outlier_margin=1,
    roof_outlier_margin=1.5,
    roof_outlier_neighbors=5,
    roof_ransac_margin=3.0,
    roof_ransac_iterations=150,
):
    builder_citymodel = builder_datamodel.create_builder_citymodel(citymodel)
    builder_pointcloud = builder_datamodel.create_builder_pointcloud(pointcloud)
    builder_citymodel = _pybuilder.extractRoofPoints(
        builder_citymodel,
        builder_pointcloud,
        ground_padding,
        ground_outlier_margin,
        roof_outlier_margin,
        roof_outlier_neighbors,
        roof_ransac_margin,
        roof_ransac_iterations,
    )
    # start_time = time()
    for citymodel_building, builder_buildings in zip(
        citymodel.buildings, builder_citymodel.buildings
    ):
        citymodel_building.roofpoints.points = np.array(
            [[p.x, p.y, p.z] for p in builder_buildings.roof_points]
        )
        ground_points = np.array(
            [[p.x, p.y, p.z] for p in builder_buildings.ground_points]
        )
        if len(ground_points) > 0:
            ground_z = ground_points[:, 2]
            citymodel_building.ground_level = np.percentile(ground_z, 50)

    # print("conver cm time: ", time() - start_time)
    return citymodel


def calculate_building_heights(
    citymodel, roof_percentile=0.9, min_height=2.5, overwrite=False
):
    for building in citymodel.buildings:
        if building.height > 0 and not overwrite:
            continue
        if len(building.roofpoints) == 0:
            building.height = min_height
            continue
        z_values = building.roofpoints.points[:, 2]
        roof_top = np.percentile(z_values, roof_percentile * 100)

        if building.ground_level == 0:
            if len(citymodel.terrain.shape) == 2:
                footprint_center = building.footprint.centroid
                ground_height = citymodel.terrain.get_value(
                    footprint_center.x, footprint_center.y
                )
                building.ground_level = ground_height
        height = roof_top - building.ground_level
        if height < min_height:
            height = min_height
        building.height = height

    return citymodel


def build_mesh(
    city_model: model.CityModel,
    mesh_resolution,
    domain_height,
    min_building_dist,
    min_vertex_dist,
    debug=False,
) -> Tuple[model.VolumeMesh, model.Mesh]:
    builder_cm = builder_datamodel.create_builder_citymodel(city_model)
    bounds = (
        city_model.bounds.xmin,
        city_model.bounds.ymin,
        city_model.bounds.xmax,
        city_model.bounds.ymax,
    )
    print(f"Building mesh with bounds {bounds}")
    print(f"dem bounds {city_model.terrain.bounds}")
    simple_cm = _pybuilder.SimplifyCityModel(
        builder_cm, bounds, min_building_dist, min_vertex_dist
    )
    # simple_cm = _pybuilder.CleanCityModel(simple_cm, min_vertex_dist)

    builder_dem = builder_datamodel.raster_to_builder_gridfield(city_model.terrain)

    # Step 3.1: Generate 2D mesh
    mesh_2D = _pybuilder.GenerateMesh2D(
        simple_cm,
        bounds,
        mesh_resolution,
    )

    if debug:
        builder_datamodel.builder_mesh2D_to_mesh(mesh_2D).save("mesh_step3.1.vtu")

    # Step 3.2: Generate 3D mesh (layer 3D mesh)
    mesh_3D = _pybuilder.GenerateMesh3D(mesh_2D, domain_height, mesh_resolution)
    if debug:
        builder_datamodel.builder_mesh3D_to_volume_mesh(mesh_3D).save(
            "meshing_step3.2.vtu"
        )

    # FIXME: Make parameters
    max_iterations = 1000
    relative_tolerance = 1e-3

    # Step 3.3: Smooth 3D mesh (set ground height)
    top_height = domain_height + city_model.terrain.data.mean()
    mesh_3D = _pybuilder.smooth_volume_mesh(
        mesh_3D,
        simple_cm,
        builder_dem,
        top_height,
        False,
        max_iterations,
        relative_tolerance,
    )

    if debug:
        builder_datamodel.builder_mesh3D_to_volume_mesh(mesh_3D).save(
            "meshing_step3.3.vtu"
        )

    # Step 3.4: Trim 3D mesh (remove building interiors)
    mesh_3D = _pybuilder.TrimMesh3D(mesh_3D, mesh_2D, simple_cm)
    if debug:
        builder_datamodel.builder_mesh3D_to_volume_mesh(mesh_3D).save(
            "meshing_step3.4.vtu"
        )

    # Step 3.5: Smooth 3D mesh (set ground and building heights)
    mesh_3D = _pybuilder.smooth_volume_mesh(
        mesh_3D,
        simple_cm,
        builder_dem,
        top_height,
        True,
        max_iterations,
        relative_tolerance,
    )

    if debug:
        builder_datamodel.builder_mesh3D_to_volume_mesh(mesh_3D).save(
            "meshing_step3.5.vtu"
        )
    surface_3d = _pybuilder.ExtractBoundary3D(mesh_3D)

    # Step 3.6: Convert back to dtcc format
    dtcc_volume = builder_datamodel.builder_mesh3D_to_volume_mesh(mesh_3D)
    dtcc_surface = builder_datamodel.builder_surface_mesh_to_mesh(surface_3d)

    return dtcc_volume, dtcc_surface


def build_surface_meshes(
    city_model: model.CityModel, min_building_dist, min_vertex_dist, mesh_resolution
):
    bounds = (
        city_model.bounds.xmin,
        city_model.bounds.ymin,
        city_model.bounds.xmax,
        city_model.bounds.ymax,
    )
    builder_cm = builder_datamodel.create_builder_citymodel(city_model)
    simple_cm = _pybuilder.SimplifyCityModel(
        builder_cm, bounds, min_building_dist, min_vertex_dist
    )
    # simple_cm = _pybuilder.CleanCityModel(simple_cm, min_vertex_dist)

    builder_dem = builder_datamodel.raster_to_builder_gridfield(city_model.terrain)

    surfaces = _pybuilder.GenerateSurface3D(simple_cm, builder_dem, mesh_resolution)

    ground_surface = surfaces[0]
    building_surfaces = surfaces[1:]
    building_surfaces = _pybuilder.MergeSurfaces3D(building_surfaces)

    dtcc_ground_surface = builder_datamodel.builder_surface_mesh_to_mesh(ground_surface)
    dtcc_building_surfaces = builder_datamodel.builder_surface_mesh_to_mesh(
        building_surfaces
    )

    return dtcc_ground_surface, dtcc_building_surfaces
