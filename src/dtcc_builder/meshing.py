from . import _dtcc_builder
from dtcc_builder.model import create_builder_polygon, builder_mesh_to_mesh
from dtcc_model import Building, City, PointCloud


def extrude_building(building: Building, resolution=5, ground_to_zero: bool = False):
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
    else:
        ground = building.ground_level
    builder_mesh = _dtcc_builder.extrude_footprint(
        builder_polygon, resolution, ground, building.height
    )
    mesh = builder_mesh_to_mesh(builder_mesh)
    return mesh
