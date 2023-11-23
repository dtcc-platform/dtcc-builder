from . import _dtcc_builder
from dtcc_builder.model import (
    create_builder_polygon,
    builder_mesh_to_mesh,
    mesh_to_builder_mesh,
    create_builder_city,
    raster_to_builder_gridfield,
)
from time import time
from dtcc_model import Building, City, PointCloud, Mesh
from .logging import info, warning, error, debug


def terrain_mesh(
    city: City, max_mesh_size, min_mesh_angle, smoothing, include_footprints=True
):
    if city.terrain is None or city.terrain.data.shape[0] == 0:
        raise ValueError("City has no terrain data. Please compute terrain first.")
    if include_footprints:
        merged_city = city.merge_buildings()
        builder_city = create_builder_city(merged_city)
    else:
        builder_city = _dtcc_builder.City()
    builder_dem = raster_to_builder_gridfield(city.terrain)

    ground_mesh = _dtcc_builder.build_terrain_mesh(
        builder_city, builder_dem, max_mesh_size, min_mesh_angle, smoothing
    )

    ground_mesh = builder_mesh_to_mesh(ground_mesh)

    return ground_mesh


def extrude_building(
    building: Building,
    max_mesh_size,
    min_mesh_angle,
    ground_to_zero: bool = False,
    cap_base: bool = False,
    per_floor=False,
):
    """
    Extrude the given building to its height.

    Parameters
    ----------
    `building` : dtcc_model.Building
        The building to extrude.
    `ground_to_zero` : bool, optional
        Whether to set the ground level to zero, by default False.

    Returns
    -------
    `dtcc_model.Building`
        The extruded building.
    """
    builder_polygon = create_builder_polygon(building.footprint)
    if ground_to_zero:
        ground = 0
        height = building.height

    else:
        ground = building.ground_level
        height = ground + building.height
    if (not per_floor) or building.floors <= 1:
        builder_mesh = _dtcc_builder.extrude_footprint(
            builder_polygon, max_mesh_size, min_mesh_angle, ground, height, cap_base
        )
        mesh = builder_mesh_to_mesh(builder_mesh)
    else:
        floor_height = building.height / building.floors
        floor_meshes = []
        for i in range(building.floors):
            if i == 0:
                cap_base = cap_base
            else:
                cap_base = False
            floor_base = ground + (i * floor_height)
            floor_mesh = _dtcc_builder.extrude_footprint(
                builder_polygon,
                max_mesh_size,
                min_mesh_angle,
                floor_base,
                floor_base + floor_height,
                cap_base,
            )
            floor_meshes.append(floor_mesh)
        merged_mesh = _dtcc_builder.merge_meshes(floor_meshes, False)
        mesh = builder_mesh_to_mesh(merged_mesh)

    return mesh


def building_meshes(
    city: City,
    max_mesh_size,
    min_mesh_angle,
    ground_to_zero: bool = False,
    cap_base: bool = False,
    per_floor=False,
    single_mesh=False,
):
    meshes = []
    for building in city.buildings:
        meshes.append(
            extrude_building(
                building,
                max_mesh_size,
                min_mesh_angle,
                ground_to_zero,
                cap_base,
                per_floor,
            )
        )
    if single_mesh:
        meshes = merge_meshes(meshes, False)
    return meshes


def city_surface_mesh(
    city: City,
    max_mesh_size,
    min_mesh_angle,
    smoothing=0,
    merge_meshes=True,
    merge_buildings=True,
    min_building_detail=1,
    min_building_area=15,
    height_merge_strategy="area_weighted",
):
    if city.terrain is None or city.terrain.data.shape[0] == 0:
        raise ValueError("City has no terrain data. Please compute terrain first.")
    start_time = time()
    if merge_buildings:
        merged_city = city.merge_buildings(
            min_building_detail,
            min_building_area,
            height_merge_strategy=height_merge_strategy,
        )  # needed for triangulation
        merged_city = merged_city.fix_building_clearance(min_building_detail, 30)
        debug(f"merge buildings took {time() - start_time} seconds")
    else:  # we've already merged the buildings (hopefully...or we crash)
        merged_city = city
    start_time = time()
    builder_city = create_builder_city(merged_city)
    debug(f"create builder city took {time() - start_time} seconds")

    builder_dem = raster_to_builder_gridfield(city.terrain)

    start_time = time()
    meshes = _dtcc_builder.build_city_surface_mesh(
        builder_city,
        builder_dem,
        max_mesh_size,
        min_mesh_angle,
        smoothing,
        merge_meshes,
    )
    debug(f"builder city surface mesh took {time() - start_time} seconds")
    start_time = time()
    if merge_meshes:
        # meshes contain only one merged mesh
        surface_mesh = builder_mesh_to_mesh(meshes[0])
    else:
        surface_mesh = [builder_mesh_to_mesh(mesh) for mesh in meshes]
    debug(f"convert builder mesh to mesh took {time() - start_time} seconds")
    return surface_mesh


def merge_meshes(meshes: [Mesh], weld=False) -> Mesh:
    builder_meshes = [mesh_to_builder_mesh(mesh) for mesh in meshes]
    merged_mesh = _dtcc_builder.merge_meshes(builder_meshes, weld)
    mesh = builder_mesh_to_mesh(merged_mesh)
    return mesh
