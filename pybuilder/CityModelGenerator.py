import _pybuilder

import fiona

def buildingBounds(shp_footprint_file, buffer=0):
    with fiona.open(shp_footprint_file) as c:
        bbox = c.bounds
    if buffer!=0:
        px,py,qx,qy = bbox
        bbox = (px-buffer, py-buffer, qx+buffer, qy+buffer)
    return bbox

def generateCityModel(shp_footprint_file, parameteres, bounds=None):
    if bounds is None:
        bounds = buildingBounds(shp_footprint_file,parameteres["DomainMargin"])
    city_model = _pybuilder.GenerateCityModel(str(shp_footprint_file), bounds,
                                              float(parameteres["MinBuildingDistance"]), 
                                              float(parameteres["MinBuildingSize"])) 
    return city_model


