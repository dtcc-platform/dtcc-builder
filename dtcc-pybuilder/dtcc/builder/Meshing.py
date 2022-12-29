from dtcc.builder import _pybuilder
import os
from typing import List, Tuple
from dtcc.builder import CityModel
from dtcc.builder import ElevationModel


def generate_mesh2D(
    city_model: CityModel,
    bbox: Tuple[float, float, float, float],
    resolution,
) -> _pybuilder.Mesh2D:
    return _pybuilder.GenerateMesh2D(city_model._builder_cm, bbox, resolution)


def generate_mesh3D(
    mesh: _pybuilder.Mesh2D, domain_height, mesh_resolution
) -> _pybuilder.Mesh3D:
    return _pybuilder.GenerateMesh3D(mesh, domain_height, mesh_resolution)


def smooth_mesh3D(
    mesh: _pybuilder.Mesh3D,
    city_model: CityModel,
    dtm: ElevationModel,
    top_height: float,
    fix_buildings: bool,
):
    return _pybuilder.SmoothMesh3D(
        mesh, city_model._builder_cm, dtm._grid_field, top_height, fix_buildings
    )


def generate_surface3D(city_model: CityModel, dtm: ElevationModel, resolution):
    return _pybuilder.GenerateSurfaces3D(
        city_model._builder_cm, dtm._grid_field, resolution
    )


def trim_mesh3D(
    mesh: _pybuilder.Mesh3D,
    mesh2D: _pybuilder.Mesh2D,
    city_model: CityModel,
    num_layers,
) -> _pybuilder.Mesh3D:
    return _pybuilder.TrimMesh3D(mesh, mesh2D, city_model._builder_cm, num_layers)


def extract_boundary3D(mesh: _pybuilder.Mesh3D) -> _pybuilder.Surface3D:
    return _pybuilder.ExtractBoundary3D(mesh)


def extract_open_surface3D(boundary: _pybuilder.Surface3D) -> _pybuilder.Surface3D:
    return _pybuilder.ExtractOpenSurface3D(boundary)


def merge_surfaces3D(surfaces: List[_pybuilder.Surface3D]) -> _pybuilder.Surface3D:
    return _pybuilder.MergeSurfaces3D(surfaces)


def write_VTK_mesh3D(mesh: _pybuilder.Mesh3D, out_file):
    return _pybuilder.WriteVTKMesh3D(mesh, str(out_file))


def write_VTK_mesh2D(mesh: _pybuilder.Mesh2D, out_file):
    return _pybuilder.WriteVTKMesh2D(mesh, str(out_file))


def write_VTK_surface3D(surface: _pybuilder.Surface3D, out_file):
    return _pybuilder.WriteVTKSurface3D


def write_Protobuf_surface3D(surface: _pybuilder.Surface3D, out_file):
    with open(out_file, "wb") as f:
        pb = _pybuilder.convertSurface3DToProtobuf(surface)
        f.write(pb)


def write_surface(surface, file_name, format="", y_up=None):
    supported_formats = ["obj", "stl", "gltf"]
    y_up_format = ["obj", "gltf"]
    if not format:
        format = os.path.splitext(file_name)[-1][1:].lower()
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
    return _pybuilder.WriteSurface3D(surface, str(file_name), format, y_up)
