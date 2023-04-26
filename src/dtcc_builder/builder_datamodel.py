from dtcc_builder import _pybuilder
import dtcc_model as model
from dtcc_model import CityModel, PointCloud, Raster
from typing import Union
import numpy as np


def create_builder_pointcloud(
    pc: Union[PointCloud, np.ndarray]
) -> _pybuilder.PointCloud:
    if isinstance(pc, np.ndarray):
        return _pybuilder.createBuilderPointCloud(
            pc, np.empty(0), np.empty(0), np.empty(0)
        )
    else:
        return _pybuilder.createBuilderPointCloud(
            pc.points, pc.classification, pc.return_number, pc.number_of_returns
        )


def create_builder_citymodel(cm: CityModel):
    building_shells = [
        list(building.footprint.exterior.coords[:-1]) for building in cm.buildings
    ]

    uuids = [building.uuid for building in cm.buildings]
    heights = [building.height for building in cm.buildings]
    ground_levels = [building.ground_level for building in cm.buildings]
    origin = cm.origin
    return _pybuilder.createBuilderCityModel(
        building_shells, uuids, heights, ground_levels, origin
    )


def raster_to_builder_gridfield(raster: Raster):
    return _pybuilder.createBuilderGridField(
        raster.data.flatten(),
        raster.bounds.tuple,
        raster.width,
        raster.height,
        abs(raster.cell_size[0]),
        abs(raster.cell_size[1]),
    )


def builder_mesh2D_to_mesh(mesh2D: _pybuilder.Mesh2D):
    mesh = model.Mesh()
    mesh.vertices = np.array([[v.x, v.y, 0] for v in mesh2D.Vertices])
    mesh.faces = np.array([[f.v0, f.v1, f.v2] for f in mesh2D.Cells])
    mesh.normals = np.array([[0, 0, 1] for n in range(len(mesh2D.Cells))])
    return mesh


def builder_mesh3D_to_volume_mesh(mesh3D: _pybuilder.Mesh3D):
    volume_mesh = model.VolumeMesh()
    volume_mesh.vertices = np.array([[v.x, v.y, v.z] for v in mesh3D.Vertices])
    volume_mesh.cells = np.array([[c.v0, c.v1, c.v2, c.v3] for c in mesh3D.Cells])
    volume_mesh.markers = np.array(mesh3D.Markers)
    return volume_mesh


def builder_surface_mesh_to_mesh(surfaceMesh: _pybuilder.Surface3D):
    surface_mesh = model.Mesh()
    surface_mesh.vertices = np.array([[v.x, v.y, v.z] for v in surfaceMesh.Vertices])
    surface_mesh.faces = np.array([[f.v0, f.v1, f.v2] for f in surfaceMesh.Faces])
    surface_mesh.normals = np.array([[n.x, n.y, n.z] for n in surfaceMesh.Normals])
    return surface_mesh
