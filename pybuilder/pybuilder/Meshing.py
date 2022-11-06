import _pybuilder
import os
import numpy
from typing import List, Tuple


def generateMesh2D(
    city_model: _pybuilder.CityModel,
    bbox: Tuple[float, float, float, float],
    resolution,
) -> _pybuilder.Mesh2D:
    return _pybuilder.GenerateMesh2D(city_model, bbox, resolution)


def generateMesh3D(
    mesh: _pybuilder.Mesh2D, domain_height, mesh_resolution
) -> _pybuilder.Mesh3D:
    return _pybuilder.GenerateMesh3D(mesh, domain_height, mesh_resolution)


def smoothMesh3D(
    mesh: _pybuilder.Mesh3D,
    city_model: _pybuilder.CityModel,
    dem: _pybuilder.GridField2D,
    top_height: float,
    fix_buildings: bool,
):
    return _pybuilder.SmoothMesh3D(mesh, city_model, dem, top_height, fix_buildings)


def generateSurface3D(
    city_model: _pybuilder.CityModel, dtm: _pybuilder.GridField2D, resolution
):
    return _pybuilder.GenerateSurfaces3D(city_model, dtm, resolution)


def trimMesh3D(
    mesh: _pybuilder.Mesh3D,
    mesh2D: _pybuilder.Mesh2D,
    city_model: _pybuilder.CityModel,
    num_layers,
) -> _pybuilder.Mesh3D:
    return _pybuilder.TrimMesh3D(mesh, mesh2D, city_model, num_layers)


def extractBoundary3D(mesh: _pybuilder.Mesh3D) -> _pybuilder.Surface3D:
    return _pybuilder.ExtractBoundary3D(mesh)


def extractOpenSurface3D(boundary: _pybuilder.Surface3D) -> _pybuilder.Surface3D:
    return _pybuilder.ExtractOpenSurface3D(boundary)


def mergeSurface3D(surfaces: List[_pybuilder.Surface3D]) -> _pybuilder.Surface3D:
    return _pybuilder.MergeSurface3D(surfaces)
