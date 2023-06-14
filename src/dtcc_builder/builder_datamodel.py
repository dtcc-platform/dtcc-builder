from dtcc_builder import _dtcc_builder
import dtcc_model as model
from dtcc_model import CityModel, PointCloud, Raster
from typing import Union
import numpy as np


def create_builder_pointcloud(
    pc: Union[PointCloud, np.ndarray]
) -> _dtcc_builder.PointCloud:
    if isinstance(pc, np.ndarray):
        return _dtcc_builder.createBuilderPointCloud(
            pc, np.empty(0), np.empty(0), np.empty(0)
        )
    else:
        return _dtcc_builder.createBuilderPointCloud(
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
    return _dtcc_builder.createBuilderCityModel(
        building_shells, uuids, heights, ground_levels, origin
    )


def raster_to_builder_gridfield(raster: Raster):
    return _dtcc_builder.createBuilderGridField(
        raster.data.flatten(),
        raster.bounds.tuple,
        raster.width,
        raster.height,
        abs(raster.cell_size[0]),
        abs(raster.cell_size[1]),
    )


def builder_mesh_to_mesh(_mesh: _dtcc_builder.Mesh):
    mesh = model.Mesh()
    mesh.vertices = np.array([[v.x, v.y, v.z] for v in _mesh.Vertices])
    mesh.faces = np.array([[f.v0, f.v1, f.v2] for f in _mesh.Faces])
    mesh.normals = np.array([[n.x, n.y, n.z] for n in _mesh.Normals])
    return mesh


def builder_volume_mesh_to_volume_mesh(_volume_mesh: _dtcc_builder.Mesh):
    volume_mesh = model.VolumeMesh()
    volume_mesh.vertices = np.array([[v.x, v.y, v.z] for v in _volume_mesh.Vertices])
    volume_mesh.cells = np.array([[c.v0, c.v1, c.v2, c.v3] for c in _volume_mesh.Cells])
    volume_mesh.markers = np.array(_volume_mesh.Markers)
    return volume_mesh


