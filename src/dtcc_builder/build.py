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
    ground_point_idx = np.where(np.isin(pointcloud.classification, [2, 9]))[0]
    ground_points = pointcloud.points[ground_point_idx]
    if isinstance(bounds, model.Bounds):
        bounds = bounds.tuple
    dem = points2grid(ground_points, cell_size, bounds, window_size=window_size)

    dem_raster = model.Raster()
    dem_raster.data = dem
    dem_raster.nodata = 0
    dem_raster.georef = Affine.translation(*bounds[:2]) * Affine.scale(
        cell_size, cell_size
    )
    io.process.raster.fill_holes(dem_raster)
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
    start_time = time()
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

    print("conver cm time: ", time() - start_time)
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
    num_layers,
    min_building_dist,
    min_vertex_dist,
) -> Tuple[model.VolumeMesh, model.Mesh]:
    builder_cm = builder_datamodel.create_builder_citymodel(city_model)
    bounds = (
        city_model.bounds.xmin,
        city_model.bounds.ymin,
        city_model.bounds.xmax,
        city_model.bounds.ymax,
    )
    simple_cm = _pybuilder.SimplifyCityModel(
        builder_cm, bounds, min_building_dist, min_vertex_dist
    )
    simple_cm = _pybuilder.CleanCityModel(simple_cm, min_vertex_dist)

    builder_dem = builder_datamodel.raster_to_builder_gridfield(city_model.terrain)

    mesh_2D = _pybuilder.GenerateMesh2D(
        simple_cm,
        bounds,
        mesh_resolution,
    )
    mesh_3D = _pybuilder.GenerateMesh3D(mesh_2D, domain_height, mesh_resolution)

    mesh_3D = _pybuilder.SmoothMesh3D(mesh_3D, simple_cm, builder_dem, False)

    mesh_3D = _pybuilder.TrimMesh3D(mesh_3D, mesh_2D, simple_cm, num_layers)

    mesh_3D = _pybuilder.SmoothMesh3D(mesh_3D, simple_cm, builder_dem, True)

    surface_3d = _pybuilder.ExtractBoundary3D(mesh_3D)


def build_surface_meshes(
    building_footprint_path, pointcloud_path, parameters=None, return_protobuf=True
):
    """Create a surface mesh from a directory of point clouds and a shapefile of building footprints.
    Args:
    @param building_footprint_path: Path to a shapefile containing building footprints
    @param pointcloud_path: Path to a pointcloud file or a directory containing point clouds
    @param parameters: parameters dict or Path to a JSON file containing parameters for the city model.
    If None default parameters will be used.
    @param return_protobuf: If True, the city model will be returned as a protobuf string

    @return: A tuple containing the ground surface and building surfaces as a protobuf string
    (if return_protobuf is True) or as builder objects (if return_protobuf is False)
    """
    if type(parameters) is dict:
        p = parameters
    else:
        p = builder_datamodel.load_parameters(parameters)

    cm, dtm = build_citymodel(building_footprint_path, pointcloud_path, p, False)
    if not cm.simplified:
        cm.simplify_citymodel(dtm.bounds)
    surfaces = builder_datamodel.Meshing.generate_surface3D(
        cm, dtm, p["MeshResolution"]
    )
    ground_surface = surfaces[0]
    building_surfaces = surfaces[1:]
    building_surfaces = builder_datamodel.Meshing.merge_surfaces3D(building_surfaces)
    if return_protobuf:
        gs = model.Mesh()
        bs = model.Mesh()
        gs.ParseFromString(_pybuilder.convertSurface3DToProtobuf(ground_surface))
        bs.ParseFromString(_pybuilder.convertSurface3DToProtobuf(building_surfaces))
        ground_surface = gs
        building_surfaces = bs
    return ground_surface, building_surfaces


def build_volume_mesh(
    building_footprint_path, pointcloud_path, parameters=None, return_protobuf=True
):
    """Create a volume mesh from a directory of point clouds and a shapefile of building footprints.
    Args:
    @param building_footprint_path: Path to a shapefile containing building footprints
    @param pointcloud_path: Path to a pointcloud file or a directory containing point clouds
    @param parameters: parameters dict or Path to a JSON file containing parameters for the city model.
    If None default parameters will be used.
    @param return_protobuf: If True, the city model will be returned as a protobuf string

    @return: volume mesh and surface mesh of the city model as a protobuf string
    (if return_protobuf is True) or as builder objects (if return_protobuf is False)
    """
    if type(parameters) is dict:
        p = parameters
    else:
        p = builder_datamodel.load_parameters(parameters)

    cm, dtm = build_citymodel(building_footprint_path, pointcloud_path, p, False)
    if not cm.simplified:
        cm.simplify_citymodel(dtm.bounds)

    # Step 3.1: Generate 2D mesh
    mesh_2D = builder_datamodel.Meshing.generate_mesh2D(
        cm, dtm.bounds, p["MeshResolution"]
    )
    if p["Debug"] and p["WriteVTK"]:
        pass
        # builder.Meshing.write_surface(
        #     mesh_2D, p["OutputDirectory"] / "Step31Mesh.vtu")

    # Step 3.2: Generate 3D mesh (layer 3D mesh)
    mesh_3D = builder_datamodel.Meshing.generate_mesh3D(
        mesh_2D, p["DomainHeight"], p["MeshResolution"]
    )
    if p["Debug"] and p["WriteVTK"]:
        mesh_3D_boundary = builder_datamodel.Meshing.extract_boundary3D(mesh_3D)
        builder_datamodel.Meshing.write_mesh(
            mesh_3D_boundary, p["OutputDirectory"] / "Step32Boundary.vtu"
        )
        builder_datamodel.Meshing.write_mesh(
            mesh_3D, p["OutputDirectory"] / "Step32Mesh.vtu"
        )

    # Step 3.3: Smooth 3D mesh (set ground height)
    top_height = dtm.mean() + p["DomainHeight"]
    mesh_3D = builder_datamodel.Meshing.smooth_mesh3D(
        mesh_3D, cm, dtm, top_height, False
    )

    if p["Debug"] and p["WriteVTK"]:
        mesh_3D_boundary = builder_datamodel.Meshing.extract_boundary3D(mesh_3D)
        builder_datamodel.Meshing.write_mesh(
            mesh_3D_boundary, p["OutputDirectory"] / "Step33Boundary.vtu"
        )
        builder_datamodel.Meshing.write_mesh(
            mesh_3D, p["OutputDirectory"] / "Step33Mesh.vtu"
        )

    # Step 3.4: Trim 3D mesh (remove building interiors)
    mesh_3D = builder_datamodel.Meshing.trim_mesh3D(
        mesh_3D, mesh_2D, cm, mesh_3D.numLayers
    )

    # Step 3.5: Smooth 3D mesh (set ground and building heights)
    mesh_3D = builder_datamodel.Meshing.smooth_mesh3D(
        mesh_3D, cm, dtm, top_height, True
    )
    mesh_3D_boundary = builder_datamodel.Meshing.extract_boundary3D(mesh_3D)
    if p["Debug"] and p["WriteVTK"]:
        builder_datamodel.Meshing.write_mesh(
            mesh_3D_boundary, p["OutputDirectory"] / "Step35Boundary.vtu"
        )
        builder_datamodel.Meshing.write_mesh(
            mesh_3D, p["OutputDirectory"] / "Step35Mesh.vtu"
        )

    if return_protobuf:
        m3d = model.VolumeMesh()
        m3d_b = model.Mesh()
        m3d.ParseFromString(_pybuilder.convertMesh3DToProtobuf(mesh_3D))
        m3d_b.ParseFromString(_pybuilder.convertSurface3DToProtobuf(mesh_3D_boundary))
        mesh_3D = m3d
        mesh_3D_boundary = m3d_b
    return mesh_3D, mesh_3D_boundary
