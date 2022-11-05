import _pybuilder

import numpy
import fiona
from shapely.geometry import Polygon
from typing import List, Tuple


def buildingBounds(shp_footprint_file, buffer=0):
    with fiona.open(shp_footprint_file) as c:
        bbox = c.bounds
    if buffer != 0:
        px, py, qx, qy = bbox
        bbox = (px - buffer, py - buffer, qx + buffer, qy + buffer)
    return bbox


def generateCityModel(shp_footprint_file, parameteres, bounds=None):
    if bounds is None:
        bounds = buildingBounds(shp_footprint_file, parameteres["DomainMargin"])
    city_model = _pybuilder.GenerateCityModel(
        str(shp_footprint_file),
        bounds,
        float(parameteres["MinBuildingDistance"]),
        float(parameteres["MinBuildingSize"]),
    )
    return city_model


def setOrigin(
    city_model: _pybuilder.CityModel, origin: Tuple[float, float]
) -> _pybuilder.CityModel:
    return _pybuilder.SetCityModelOrigin(city_model, origin)


def cleanCityModel(
    city_model: _pybuilder.CityModel, min_vert_distance
) -> _pybuilder.CityModel:
    return _pybuilder.CleanCityModel(city_model, min_vert_distance)


def extractBuildingPoints(
    city_model: _pybuilder.CityModel,
    point_cloud: _pybuilder.PointCloud,
    ground_margins: float,
    groundOutlierMargin: float,
) -> _pybuilder.CityModel:
    return _pybuilder.ExtractBuildingPoints(
        city_model, point_cloud, ground_margins, groundOutlierMargin
    )


def buildingPointsRANSACOutlierRemover(
    city_model: _pybuilder.CityModel, outlier_margin: float, interations: int
) -> _pybuilder.CityModel:
    return _pybuilder.BuildingPointsRANSACOutlierRemover(
        city_model, outlier_margin, interations
    )


def buildingPointsStatisticalOutlierRemover(
    city_model: _pybuilder.CityModel, neighbors: int, outlier_margin: float
) -> _pybuilder.CityModel:
    return _pybuilder.BuildingPointsOutlierRemover(
        city_model, neighbors, outlier_margin
    )


def computeBuildingHeights(
    city_model: _pybuilder.CityModel,
    dtm: _pybuilder.GridField2D,
    groundPercentile: float,
    roofPercentile: float,
) -> _pybuilder.CityModel:
    return _pybuilder.ComputeBuildingHeights(
        city_model, dtm, groundPercentile, roofPercentile
    )


def getBuildingRoofPoints(building: _pybuilder.Building):
    roof_pts = building.roofPoints
    roof_pts = [[p.x, p.y, p.z] for p in roof_pts]
    return numpy.array(roof_pts)


def getBuildingFootprint(building: _pybuilder.Building):
    footprint = building.footprint


def toJSON(city_model: _pybuilder.CityModel, outfile):
    _pybuilder.WriteCityModelJSON(city_model, outfile)


def fromJSON(infile) -> _pybuilder.CityModel:
    return _pybuilder.ReadCityModelJSON(infile)
