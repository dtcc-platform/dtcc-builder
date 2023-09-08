from . import _dtcc_builder
from dtcc_builder.model import create_builder_polygon, builder_mesh_to_mesh
from dtcc_model import Building, City, PointCloud


def extrude_building(
    building: Building,
    resolution=5,
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
    if not per_floor:
        builder_mesh = _dtcc_builder.extrude_footprint(
            builder_polygon, resolution, ground, height, cap_base
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
                resolution,
                floor_base,
                floor_base + floor_height,
                cap_base,
            )
            floor_meshes.append(floor_mesh)
        merged_mesh = _dtcc_builder.merge_meshes(floor_meshes)
        mesh = builder_mesh_to_mesh(merged_mesh)

    return mesh
