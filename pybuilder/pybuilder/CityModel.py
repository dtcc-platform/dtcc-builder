import _pybuilder

import numpy
import fiona
from shapely.geometry import Polygon
from typing import List, Tuple


def building_bounds(shp_footprint_file, buffer=0):
    with fiona.open(shp_footprint_file) as c:
        bbox = c.bounds
    if buffer != 0:
        px, py, qx, qy = bbox
        bbox = (px - buffer, py - buffer, qx + buffer, qy + buffer)
    return bbox


def generate_citymodel(shp_footprint_file, parameteres, bounds=None):
    if bounds is None:
        bounds = building_bounds(shp_footprint_file, parameteres["DomainMargin"])
    city_model = _pybuilder.GenerateCityModel(
        str(shp_footprint_file),
        bounds,
        float(parameteres["MinBuildingDistance"]),
        float(parameteres["MinBuildingSize"]),
    )
    return city_model


def set_origin(
    city_model: _pybuilder.CityModel, origin: Tuple[float, float]
) -> _pybuilder.CityModel:
    return _pybuilder.SetCityModelOrigin(city_model, origin)


def clean_citymodel(
    city_model: _pybuilder.CityModel, min_vert_distance
) -> _pybuilder.CityModel:
    return _pybuilder.CleanCityModel(city_model, min_vert_distance)


def extract_building_points(
    city_model: _pybuilder.CityModel,
    point_cloud: _pybuilder.PointCloud,
    ground_margins: float,
    groundOutlierMargin: float,
) -> _pybuilder.CityModel:
    return _pybuilder.ExtractBuildingPoints(
        city_model, point_cloud, ground_margins, groundOutlierMargin
    )


def building_points_RANSAC_outlier_remover(
    city_model: _pybuilder.CityModel, outlier_margin: float, interations: int
) -> _pybuilder.CityModel:
    return _pybuilder.BuildingPointsRANSACOutlierRemover(
        city_model, outlier_margin, interations
    )


def building_points_statistical_outlier_remover(
    city_model: _pybuilder.CityModel, neighbors: int, outlier_margin: float
) -> _pybuilder.CityModel:
    return _pybuilder.BuildingPointsOutlierRemover(
        city_model, neighbors, outlier_margin
    )


def compute_building_heights(
    city_model: _pybuilder.CityModel,
    dtm: _pybuilder.GridField2D,
    groundPercentile: float,
    roofPercentile: float,
) -> _pybuilder.CityModel:
    return _pybuilder.ComputeBuildingHeights(
        city_model, dtm, groundPercentile, roofPercentile
    )


def get_building_roof_points(building: _pybuilder.Building):
    roof_pts = building.roofPoints
    roof_pts = [[p.x, p.y, p.z] for p in roof_pts]
    return numpy.array(roof_pts)


def get_building_footprint(building: _pybuilder.Building):
    footprint = building.footprint


def to_JSON(city_model: _pybuilder.CityModel, outfile):
    _pybuilder.WriteCityModelJSON(city_model, outfile)


def from_JSON(infile) -> _pybuilder.CityModel:
    return _pybuilder.ReadCityModelJSON(infile)


def load_protobuf(protobuf_string: str) -> _pybuilder.CityModel:
    return _pybuilder.loadCityModelProtobuf(protobuf_string)


def to_protobuf(cm: _pybuilder.CityModel) -> str:
    return _pybuilder.convertCityModelToProtobuf(cm)
