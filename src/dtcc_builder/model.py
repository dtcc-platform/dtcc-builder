from . import _dtcc_builder
import dtcc_model as model
from shapely.geometry import Polygon
from dtcc_model import City, PointCloud, Raster
from typing import Union
import numpy as np


def create_builder_polygon(polygon: Polygon) -> _dtcc_builder.Polygon:
    """
    Create a DTCC builder Polygon object through the pybind exposed C++
    `DTCC_BUILDER::create_polygon()` function.

    Parameters
    ----------
    polygon : model.Polygon
        The input Polygon object.

    Returns
    -------
    _dtcc_builder.Polygon
        A `DTCC_BUILDER` Polygon object.

    """
    shell = list(polygon.exterior.coords[:-1])
    holes = [list(hole.coords[:-1]) for hole in polygon.interiors]

    return _dtcc_builder.create_polygon(shell, holes)


def create_builder_pointcloud(
    pc: Union[PointCloud, np.ndarray]
) -> _dtcc_builder.PointCloud:
    """
    Create a DTCC builder PointCloud object through the pybind exposed C++
    `DTCC_BUILDER::create_pointcloud()` function.

    Parameters
    ----------
    pc : Union[PointCloud, np.ndarray]
        The input PointCloud data or numpy array.

    Returns
    -------
    _dtcc_builder.PointCloud
        A `DTCC_BUILDER` PointCloud object.

    """
    if isinstance(pc, np.ndarray):
        return _dtcc_builder.create_pointcloud(
            pc, np.empty(0), np.empty(0), np.empty(0)
        )
    else:
        return _dtcc_builder.create_pointcloud(
            pc.points, pc.classification, pc.return_number, pc.num_returns
        )


def create_builder_surface(surface: model.Surface):
    """
    Create a DTCC builder Surface object through the pybind exposed C++
    `DTCC_BUILDER::create_surface()` function.

    Parameters
    ----------
    surface : model.Surface
        The input Surface object.

    Returns
    -------
    _dtcc_builder.Surface
        A `DTCC_BUILDER` Surface object.

    """
    return _dtcc_builder.create_surface(surface.vertices, surface.holes)


def create_builder_multisurface(multisurface: model.MultiSurface):
    """
    Create a DTCC builder MultiSurface object through the pybind exposed C++
    `DTCC_BUILDER::create_multisurface()` function.

    Parameters
    ----------
    multisurface : model.MultiSurface
        The input MultiSurface object.

    Returns
    -------
    _dtcc_builder.MultiSurface
        A `DTCC_BUILDER` MultiSurface object.

    """
    surfaces = [
        _dtcc_builder.create_surface(surface.vertices, surface.holes)
        for surface in multisurface.surfaces
    ]
    return _dtcc_builder.create_multisurface(surfaces)


def raster_to_builder_gridfield(raster: Raster):
    """
    Convert Raster to a DTCC builder GridField object through the pybind exposed C++
    `DTCC_BUILDER::create_gridfield()` function.

    Parameters
    ----------
    raster : Raster
        The input Raster object.

    Returns
    -------
    _dtcc_builder.GridField
        A `DTCC_BUILDER` GridField object.

    """
    # rasters start in top left corner, gridfields in bottom left
    # flip the data to match
    data = np.flipud(raster.data).flatten()
    return _dtcc_builder.create_gridfield(
        data,
        raster.bounds.tuple,
        raster.width,
        raster.height,
        abs(raster.cell_size[0]),
        abs(raster.cell_size[1]),
    )


def create_builder_city(city: City):
    """
    Create a DTCC builder City object through the pybind exposed C++

    `DTCC_BUILDER::create_city()` function.

    Parameters
    ----------
    city : City
        The input City object.

    Returns
    -------
    _dtcc_builder.City
        A `DTCC_BUILDER` City object.

    """
    building_shells = [
        list(building.footprint.exterior.coords[:-1]) for building in city.buildings
    ]

    building_holes = []
    for building in city.buildings:
        holes = [list(hole.coords[:-1]) for hole in building.footprint.interiors]
        building_holes.append(holes)

    uuids = [building.uuid for building in city.buildings]
    heights = [building.height for building in city.buildings]
    ground_levels = [building.ground_level for building in city.buildings]
    origin = city.origin
    return _dtcc_builder.create_city(
        building_shells, building_holes, uuids, heights, ground_levels, origin
    )


def mesh_to_builder_mesh(mesh: model.Mesh):
    """
    Convert a dtcc_model Mesh to a DTCC builder Mesh through the pybind exposed C++
    `DTCC_BUILDER::create_mesh()` function.

    Parameters
    ----------
    mesh : model.Mesh
        The input dtcc_model Mesh object.

    Returns
    -------
    _dtcc_builder.Mesh
        A DTCC builder Mesh object.

    """
    return _dtcc_builder.create_mesh(mesh.vertices, mesh.faces, mesh.markers)


def builder_mesh_to_mesh(_mesh: _dtcc_builder.Mesh):
    """
    Convert a DTCC builder Mesh to a dtcc_model Mesh.

    Parameters
    ----------
    _mesh : _dtcc_builder.Mesh
        The input DTCC builder Mesh object.

    Returns
    -------
    model.Mesh
        A dtcc_model Mesh object.

    """
    mesh = model.Mesh()
    vertices, faces, markers = _dtcc_builder.mesh_as_arrays(_mesh)
    mesh.vertices = vertices.reshape((-1, 3))
    mesh.faces = faces.reshape((-1, 3))
    mesh.markers = markers
    return mesh


def builder_volume_mesh_to_volume_mesh(_volume_mesh: _dtcc_builder.Mesh):
    """
    Convert a DTCC builder VolumeMesh to a dtcc_model VolumeMesh.

    Parameters
    ----------
    _volume_mesh : _dtcc_builder.Mesh
        The input DTCC builder VolumeMesh object.

    Returns
    -------
    model.VolumeMesh
        A dtcc_model VolumeMesh object.

    """
    volume_mesh = model.VolumeMesh()
    volume_mesh.vertices = np.array([[v.x, v.y, v.z] for v in _volume_mesh.vertices])
    volume_mesh.cells = np.array([[c.v0, c.v1, c.v2, c.v3] for c in _volume_mesh.cells])
    volume_mesh.markers = np.array(_volume_mesh.markers)
    return volume_mesh
