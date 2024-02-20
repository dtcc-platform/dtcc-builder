from dtcc_model import Building, GeometryType
from dtcc_model.geometry import Surface, MultiSurface
from logging import debug, info, warning, error
import numpy as np


def extrude_surface(surface: Surface, height: float) -> MultiSurface:
    """
    Extrude a surface to a given height. The height is the absolute height, not the height relative to the surface.

    Parameters
    ----------
    `surface` : dtcc_model.Surface
        The surface to extrude.
    `height` : float
        The height to extrude the surface to.

    Returns
    -------
    `dtcc_model.MultiSurface`
        The extruded surface.
    """
    extrusion = MultiSurface()
    base = surface.copy()
    base.vertices[:, 2] = height
    for hole in base.holes:
        hole[:, 2] = height
    extrusion.surfaces.append(surface)
    extrusion.surfaces.append(base)

    for i in range(surface.vertices.shape[0]):
        wall = Surface()
        wall.vertices = np.array(
            [
                surface.vertices[i],
                surface.vertices[(i + 1) % surface.vertices.shape[0]],
                base.vertices[(i + 1) % base.vertices.shape[0]],
                base.vertices[i],
            ]
        )
        extrusion.surfaces.append(wall)
    for i in range(len(surface.holes)):
        for j in range(surface.holes[i].shape[0]):
            hole = surface.holes[i]
            wall = Surface()
            wall.vertices = np.array(
                [
                    hole[j],
                    hole[(j + 1) % hole.shape[0]],
                    base.holes[i][(j + 1) % hole.shape[0]],
                    base.holes[i][j],
                ]
            )
            extrusion.surfaces.append(wall)
    return extrusion


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
