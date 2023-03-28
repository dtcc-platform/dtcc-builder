import dtcc_builder as builder
import dtcc_model as model
from dtcc_builder import _pybuilder
import dtcc_io as io
from pathlib import Path


def project_domain(p):
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


def build_citymodel(
    buildings, pc, parameters=None, bounds = (), return_protobuf=True
):
    """Create a city model from a directory of point clouds and a shapefile of building footprints.
    Args:
    @param building_footprint_path: Path to a shapefile containing building footprints
    @param pointcloud_path: Path to a pointcloud file or a directory containing point clouds
    @param parameters: parameters dict or Path to a JSON file containing parameters for the city model.
    If None default parameters will be used.
    @param return_protobuf: If True, the city model will be returned as a protobuf string

    @return: A tuple containing the city model and elevation model as a protobuf string
    (if return_protobuf is True) or as builder objects (if return_protobuf is False)
    """
    if type(parameters) is dict:
        p = parameters
    else:
        p = builder.load_parameters(parameters)



    builder_pc = builder.PointCloud(pc)
    if len(builder_pc) == 0:
        raise ValueError("No points in point cloud")
    builder_pc.set_origin(origin)
    builder_pc.remove_global_outliers(p["OutlierMargin"])
    if p["NaiveVegitationFilter"]:
        builder_pc.vegetation_filter()
    city_model = builder.CityModel(
        buildings,parameters)
    
    city_model.set_origin(origin)
    city_model.clean_citymodel()

    dtm = builder.ElevationModel(builder_pc, p["ElevationModelResolution"], [2, 9])
    # dsm = ElevationModel(pc, p["ElevationModelResolution"])
    dtm.smooth_elevation_model(p["GroundSmoothing"])

    city_model.extract_building_points(builder_pc)

    if 6 not in pc.used_classifications and p["RANSACOutlierRemover"]:
        city_model.building_points_RANSAC_outlier_remover()

    if p["StatisticalOutlierRemover"]:
        city_model.building_points_statistical_outlier_remover()

    city_model.compute_building_heights(dtm)

    if return_protobuf:
        cm = model.CityModel()
        cm_string = city_model.to_protobuf()
        cm.ParseFromString(cm_string)
        city_model = cm
        dtm = "" #dtm.to_protobuf()

    return (city_model, dtm)


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
        p = builder.load_parameters(parameters)

    cm, dtm = build_citymodel(building_footprint_path, pointcloud_path, p, False)
    if not cm.simplified:
        cm.simplify_citymodel(dtm.bounds)
    surfaces = builder.Meshing.generate_surface3D(cm, dtm, p["MeshResolution"])
    ground_surface = surfaces[0]
    building_surfaces = surfaces[1:]
    building_surfaces = builder.Meshing.merge_surfaces3D(building_surfaces)
    if return_protobuf:
        gs = model.Surface3D()
        bs = model.Surface3D()
        gs.ParseFromString(_pybuilder.convertSurface3DToProtobuf(ground_surface))
        bs.ParseFromString(_pybuilder.convertSurface3DToProtobuf(building_surfaces))
        ground_surface = gs
        building_surfaces = bs
    return ground_surface, building_surfaces


def build_mesh(
    city_model, parameters=None, return_protobuf=True
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
        p = builder.load_parameters(parameters)

    builder_cm = builder.CityModel(city_model)

    if not builder_cm.simplified:
        builder_cm.simplify_citymodel(builder_cm.bounds)

    # Step 3.1: Generate 2D mesh
    mesh_2D = builder.Meshing.generate_mesh2D(builder_cm, p["MeshResolution"])
    if p["Debug"] and p["WriteVTK"]:
        pass
        # builder.Meshing.write_surface(
        #     mesh_2D, p["OutputDirectory"] / "Step31Mesh.vtu")

    # Step 3.2: Generate 3D mesh (layer 3D mesh)
    mesh_3D = builder.Meshing.generate_mesh3D(
        mesh_2D, p["DomainHeight"], p["MeshResolution"]
    )
    if p["Debug"] and p["WriteVTK"]:
        mesh_3D_boundary = builder.Meshing.extract_boundary3D(mesh_3D)
        builder.Meshing.write_mesh(
            mesh_3D_boundary, p["OutputDirectory"] / "Step32Boundary.vtu"
        )
        builder.Meshing.write_mesh(
            mesh_3D, p["OutputDirectory"] / "Step32Mesh.vtu"
        )

    # Step 3.3: Smooth 3D mesh (set ground height)
    top_height = dtm.mean() + p["DomainHeight"]
    mesh_3D = builder.Meshing.smooth_mesh3D(mesh_3D, cm, dtm, top_height, False)

    if p["Debug"] and p["WriteVTK"]:
        mesh_3D_boundary = builder.Meshing.extract_boundary3D(mesh_3D)
        builder.Meshing.write_mesh(
            mesh_3D_boundary, p["OutputDirectory"] / "Step33Boundary.vtu"
        )
        builder.Meshing.write_mesh(
            mesh_3D, p["OutputDirectory"] / "Step33Mesh.vtu"
        )

    # Step 3.4: Trim 3D mesh (remove building interiors)
    mesh_3D = builder.Meshing.trim_mesh3D(mesh_3D, mesh_2D, cm, mesh_3D.numLayers)

    # Step 3.5: Smooth 3D mesh (set ground and building heights)
    mesh_3D = builder.Meshing.smooth_mesh3D(mesh_3D, cm, dtm, top_height, True)
    mesh_3D_boundary = builder.Meshing.extract_boundary3D(mesh_3D)
    if p["Debug"] and p["WriteVTK"]:
        builder.Meshing.write_mesh(
            mesh_3D_boundary, p["OutputDirectory"] / "Step35Boundary.vtu"
        )
        builder.Meshing.write_mesh(
            mesh_3D, p["OutputDirectory"] / "Step35Mesh.vtu"
        )

    if return_protobuf:
        m3d = model.VolumeMesh()
        m3d_b = model.Surface3D()
        m3d.ParseFromString(_pybuilder.convertVolumeMeshToProtobuf(mesh_3D))
        m3d_b.ParseFromString(_pybuilder.convertSurface3DToProtobuf(mesh_3D_boundary))
        mesh_3D = m3d
        mesh_3D_boundary = m3d_b
    return mesh_3D, mesh_3D_boundary
