from dtcc_model import Building, GeometryType
from dtcc_model.geometry import Surface, MultiSurface
from logging import debug, info, warning, error
import numpy as np

from .surface import extrude_surface


def extrude_building(
    building: Building, default_ground_height=0, height_property=""
) -> MultiSurface:
    """
    Extrudes the LOD0 representation of a building from its height to the ground leve.

    Parameters
    ----------
    `building` : dtcc_model.Building
        The building to extrude.
    `default_ground_height` : float, optional
        If building does not have a ground_height property, the default ground
        level to use, by default 0.

    Returns
    -------
    `MultiSurface`
        The extruded building.
    """
    ground_height = building.attributes.get("ground_height", default_ground_height)
    geometry = building.lod0
    if geometry is None:
        error(f"Building {building.id} has no LOD0 geometry.")
        return None
    if isinstance(geometry, Surface):
        geometry = MultiSurface(surface[geometry])
    if not isinstance(geometry, MultiSurface):
        error(f"Building {building.id} LOD0 geometry is not a MultiSurface.")
        return None

    extrusion = MultiSurface()
    for surface in geometry.surfaces:
        extrusion = extrusion.merge(extrude_surface(surface, ground_height))
    return extrusion


def build_lod1_buildings(
    buildings: [Building], default_ground_height=0, height_property="", rebuild=True
) -> [Building]:
    """
    Build the LOD1 representation of the given buildings.

    Parameters
    ----------
    `buildings` : [dtcc_model.Building]
        The buildings to build the LOD1 representation of.
    `default_ground_height` : float, optional
        If building does not have a ground_height property, the default ground
        level to use, by default 0.
    `height_property` : str, optional
        The property to use as the height of the building, by default "".
    `rebuild` : bool, optional
        Whether to rebuild the LOD1 representation if it already exists, by default True.
    Returns
    -------
    [dtcc_model.Building]
        The buildings with the LOD1 representation built.
    """
    for building in buildings:
        if building.lod1 is not None and not rebuild:
            continue
        if building.lod0 is None:
            warning(f"Building {building.id} has no LOD0 geometry.")
            continue
        geometry = extrude_building(building)
        if geometry is not None:
            building.add_geometry(geometry, GeometryType.LOD1)
        else:
            warning(f"Building {building.id} LOD1 geometry could not be built.")
    return buildings
