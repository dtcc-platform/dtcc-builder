import _pybuilder

import numpy
import fiona
from shapely.geometry import Polygon
from PointCloud import PointCloud
from ElevationModel import ElevationModel
from typing import List, Tuple
import Parameters


def building_bounds(shp_footprint_file, buffer=0):
    with fiona.open(shp_footprint_file) as c:
        bbox = c.bounds
    if buffer != 0:
        px, py, qx, qy = bbox
        bbox = (px - buffer, py - buffer, qx + buffer, qy + buffer)
    return bbox


class CityModel:
    def __init__(self, shp_footprints=None, parameters=None, bounds=None):
        if parameters is None:
            parameters = Parameters.load_parameters()
        self.parameters = parameters
        self.bounds = bounds
        self._builder_cm = None
        self.origin = (0, 0)
        self.cleaned = False
        self.extracted_points = False
        self.calculated_heights = False
        if shp_footprints is not None:
            self.generate_citymodel(shp_footprints)

    def get_buildings(self):
        if self._builder_cm is None:
            return None
        else:
            return self._builder_cm.buildings

    def set_buildings(self, b):
        return None

    buildings = property(get_buildings)

    def generate_citymodel(self, shp_footprint_file):
        if self.bounds is None:
            self.bounds = building_bounds(
                shp_footprint_file, self.parameters["DomainMargin"]
            )
        self._builder_cm = _pybuilder.GenerateCityModel(
            str(shp_footprint_file),
            self.bounds,
            float(self.parameters["MinBuildingDistance"]),
            float(self.parameters["MinBuildingSize"]),
        )

    def set_origin(self, origin: Tuple[float, float]):
        self._builder_cm = _pybuilder.SetCityModelOrigin(self._builder_cm, origin)
        self.origin = origin

    def clean_citymodel(self, min_vert_distance=None):
        if min_vert_distance is None:
            min_vert_distance = self.parameters["MinVertexDistance"]
        self._builder_cm = _pybuilder.CleanCityModel(
            self._builder_cm, min_vert_distance
        )
        self.cleaned = True

    def extract_building_points(
        self,
        point_cloud: PointCloud,
        ground_margins=None,
        ground_outlier_margin=None,
    ):

        if ground_margins is None:
            ground_margins = self.parameters["GroundMargin"]
        if ground_outlier_margin is None:
            ground_outlier_margin = self.parameters["OutlierMargin"]

        self._builder_cm = _pybuilder.ExtractBuildingPoints(
            self._builder_cm,
            point_cloud._builder_pc,
            ground_margins,
            ground_outlier_margin,
        )
        self.extracted_points = True

    def building_points_RANSAC_outlier_remover(
        self, outlier_margin=None, interations=None
    ):

        if outlier_margin is None:
            outlier_margin = self.parameters["RANSACOutlierMargin"]
        if interations is None:
            interations = self.parameters["RANSACIterations"]

        self._builder_cm = _pybuilder.BuildingPointsRANSACOutlierRemover(
            self._builder_cm, outlier_margin, interations
        )

    def building_points_statistical_outlier_remover(
        self, neighbors=None, outlier_margin=None
    ):

        if neighbors is None:
            neighbors = self.parameters["OutlierNeighbors"]
        if outlier_margin is None:
            outlier_margin = self.parameters["OutlierSTD"]
        self._builder_cm = _pybuilder.BuildingPointsOutlierRemover(
            self._builder_cm, neighbors, outlier_margin
        )

    def compute_building_heights(
        self, dtm: ElevationModel, ground_percentile=None, roof_percentile=None
    ):

        if ground_percentile is None:
            ground_percentile = self.parameters["GroundPercentile"]
        if roof_percentile is None:
            roof_percentile = self.parameters["RoofPercentile"]

        self._builder_cm = _pybuilder.ComputeBuildingHeights(
            self._builder_cm, dtm._grid_field, ground_percentile, roof_percentile
        )

    def to_JSON(self, outfile):
        _pybuilder.WriteCityModelJSON(self._builder_cm, outfile)

    def from_JSON(self, infile):
        self._builder_cm = _pybuilder.ReadCityModelJSON(infile)

    def load_protobuf(self, protobuf_string: str):
        self._builder_cm = _pybuilder.loadCityModelProtobuf(protobuf_string)

    def to_protobuf(self) -> str:
        return _pybuilder.convertCityModelToProtobuf(self._builder_cm)


###WIP
def get_building_roof_points(building: _pybuilder.Building):
    roof_pts = building.roofPoints
    roof_pts = [[p.x, p.y, p.z] for p in roof_pts]
    return numpy.array(roof_pts)


def get_building_footprint(building: _pybuilder.Building):
    footprint = building.footprint
