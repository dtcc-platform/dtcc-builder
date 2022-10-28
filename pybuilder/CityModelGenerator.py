import _pybuilder


def GenerateCityModel(shp_footprint_file, parameteres, bounds=None):
    city_model = _pybuilder.GenerateCityModel(str(shp_footprint_file), float(parameteres["MinBuildingDistance"]), float(parameteres["MinBuildingSize"])) 
    return city_model

