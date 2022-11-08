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


def mergeSurfaces3D(surfaces: List[_pybuilder.Surface3D]) -> _pybuilder.Surface3D:
    return _pybuilder.MergeSurfaces3D(surfaces)


def writeVTKMesh3D(mesh: _pybuilder.Mesh3D, out_file):
    return _pybuilder.WriteVTKMesh3D(mesh, out_file)


def writeVTKMesh2D(mesh: _pybuilder.Mesh2D, out_file):
    return _pybuilder.WriteVTKMesh2D(mesh, out_file)


def writeSurface(surface, file_name, format="obj", y_up=None):
    supported_formats = ["obj", "stl", "gltf"]
    y_up_format = ["obj", "gltf"]
    if format not in supported_formats:
        print(f"format {format} not supported, please use one of {supported_formats}")
        return False
    if y_up is None:
        if format in y_up_format:
            y_up = True
        else:
            y_up = False
    if format == "gltf":
        format = "gltf2"
    return _pybuilder.WriteSurface3D(surface, file_name, format, y_up)
