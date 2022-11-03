import _pybuilder

import numpy
import fiona


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


def cleanCityModel(city_model, min_vert_distance):
    return _pybuilder.CleanCityModel(city_model, min_vert_distance)


def extractBuildingPoints(
    city_model: _pybuilder.CityModel,
    point_cloud: _pybuilder.PointCloud,
    ground_margins,
    groundOutlierMargin,
):
    return _pybuilder.ExtractBuildingPoints(
        city_model, point_cloud, ground_margins, groundOutlierMargin
    )


def getBuildingRoofPoints(building: _pybuilder.Building):
    roof_pts = building.roofPoints
    roof_pts = [[p.x, p.y, p.z] for p in roof_pts]

    return numpy.array(roof_pts)
