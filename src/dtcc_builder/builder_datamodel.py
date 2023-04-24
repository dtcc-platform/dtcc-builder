from dtcc_builder import _pybuilder
from dtcc_model import CityModel, PointCloud
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
