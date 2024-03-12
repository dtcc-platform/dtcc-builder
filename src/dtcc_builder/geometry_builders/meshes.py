from dtcc_model import (
    Mesh,
    VolumeMesh,
    Building,
    Terrain,
    City,
    Surface,
    MultiSurface,
    GeometryType,
)

from dtcc_builder.model import (
    create_builder_polygon,
    create_builder_surface,
    create_builder_multisurface,
    builder_mesh_to_mesh,
    mesh_to_builder_mesh,
    create_builder_city,
    raster_to_builder_gridfield,
)

from dtcc_builder import _dtcc_builder

from dtcc_builder.polygons.polygons import (
    polygon_merger,
    simplify_polygon,
    remove_slivers,
    fix_clearance,
)

from dtcc_builder.building.modify import merge_building_footprints


def build_surface_mesh(
    city: City,
    lod: GeometryType = GeometryType.LOD1,
    min_building_detail: float = 1.0,
    min_building_area: float = 15.0,
    merge_buildings: bool = True,
    merge_tolerance: float = 0.5,
    max_mesh_size: float = 10.0,
    min_mesh_angle: float = 25.0,
    merge_meshes: bool = True,
    smoothing: int = 0,
) -> Mesh:
    """
    Build a mesh from the surfaces of the buildings in the city.

    Parameters
    ----------
    `city` : dtcc_model.City
        The city to build the mesh from.
    `min_building_detail` : float, optional
        The minimum detail of the buildin to resolve, by default 1.0.
    `min_building_area` : float, optional
        The smallest building to include, by default 15.0.
    `merge_buildings` : bool, optional
        merge building footprints, by default True.
    `max_mesh_size` : float, optional
        The maximum size of the mesh, by default 1.0.
    `min_mesh_angle` : float, optional
        The minimum angle of the mesh, by default 30.0.
    `merge_meshes` : bool, optional
        Whether to merge the meshes to a single mesh, by default True.

    `smoothing` : float, optional
        The smoothing of the mesh, by default 0.0.

    Returns
    -------
    `dtcc_model.Mesh`
    """
    buildings = city.buildings
    if merge_buildings:
        merged_buildings = merge_building_footprints(buildings, lod)
        building_footprints = [
            b.get_footprint(GeometryType.LOD0) for b in merged_buildings
        ]
    else:
        building_footprints = [b.get_footprint(lod) for b in buildings]

    terrain = city.terrain
    terrain_raster = terrain.raster
    if terrain_raster is None:
        ValueError("City has no terrain raster data. Please compute terrain first.")

    builder_dem = raster_to_builder_gridfield(terrain_raster)

    builder_surfaces = [
        create_builder_surface(p) for p in building_footprints if p is not None
    ]
    builder_mesh = _dtcc_builder.build_city_surface_mesh(
        builder_surfaces[:10],
        builder_dem,
        max_mesh_size,
        min_mesh_angle,
        smoothing,
        merge_meshes,
    )
    if merge_meshes:
        result_mesh = builder_mesh_to_mesh(builder_mesh[0])
    else:
        result_mesh = [builder_mesh_to_mesh(bm) for bm in builder_mesh]
    return result_mesh
