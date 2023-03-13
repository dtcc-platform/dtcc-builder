from dtcc_builder import _pybuilder
import os
from typing import List, Tuple
from dtcc_builder import CityModel
from dtcc_builder import ElevationModel
import dtcc_io as io
import dtcc_model as model


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


def write_Protobuf_surface3D(surface: _pybuilder.Surface3D, out_file):
    with open(out_file, "wb") as f:
        pb = _pybuilder.convertSurface3DToProtobuf(surface)
        f.write(pb)


def write_surface(surface: _pybuilder.Surface3D, file_name):
    pb = _pybuilder.convertSurface3DToProtobuf(surface)
    io.mesh.write(file_name, pb)


def write_volume_mesh(mesh: _pybuilder.Mesh3D, file_name):
    pb = _pybuilder.convertMesh3DToProtobuf(mesh)
    io.mesh.write(file_name, pb, volume_mesh=True)
