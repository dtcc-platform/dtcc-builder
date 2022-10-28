import _pybuilder


def GenerateCityModel(shp_footprints, bounds, minBuildingDistance, minBuildingSize):
    city_model = _pybuilder.GenerateCityModel(shp_footprints, minBuildingDistance, minBuildingSize)  



