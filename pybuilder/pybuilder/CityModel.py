import _pybuilder

import numpy
import fiona
from shapely.geometry import Polygon
from PointCloud import PointCloud
from typing import List, Tuple


def building_bounds(shp_footprint_file, buffer=0):
    with fiona.open(shp_footprint_file) as c:
        bbox = c.bounds
    if buffer != 0:
        px, py, qx, qy = bbox
        bbox = (px - buffer, py - buffer, qx + buffer, qy + buffer)
    return bbox


class CityModel:
    def __init__(self, shp_footprints, parameters, bounds):
        self.paramaters = parameters
        self.bounds = bounds
        self._builder_cm = None
        self.generate_citymodel(shp_footprints)
        self.origin = (0, 0)
        self.cleaned = False
        self.extracted_points = False
        self.calculated_heights = False

    def generate_citymodel(self, shp_footprint_file):
        if self.bounds is None:
            self.bounds = building_bounds(
                shp_footprint_file, self.parameteres["DomainMargin"]
            )
        self._builder_cm = _pybuilder.GenerateCityModel(
            str(shp_footprint_file),
            self.bounds,
            float(self.arameteres["MinBuildingDistance"]),
            float(self.parameteres["MinBuildingSize"]),
        )

    def set_origin(self, origin: Tuple[float, float]):
        self._builder_cm = _pybuilder.SetCityModelOrigin(self._builder_cm, origin)
        self.origin = origin

    def clean_citymodel(self, min_vert_distance=None):
        if min_vert_distance is None:
            min_vert_distance = self.paramaters["MinVertexDistance"]
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
            ground_margins = self.paramaters["GroundMargin"]
        if ground_outlier_margin is None:
            ground_outlier_margin = self.paramaters["OutlierMargin"]

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
            outlier_margin = self.paramaters["RANSACOutlierMargin"]
        if interations is None:
            interations = self.paramaters["RANSACIterations"]

        self._builder_cm = _pybuilder.BuildingPointsRANSACOutlierRemover(
            self._builder_cm, outlier_margin, interations
        )

    def building_points_statistical_outlier_remover(
        self, neighbors=None, outlier_margin=None
    ):

        if neighbors is None:
            neighbors = self.paramaters["OutlierNeighbors"]
        if outlier_margin is None:
            outlier_margin = self.paramaters["OutlierSTD"]
        self._builder_cm = _pybuilder.BuildingPointsOutlierRemover(
            self._builder_cm, neighbors, outlier_margin
        )

    def compute_building_heights(
        self, dtm: _pybuilder.GridField2D, ground_percentile=None, roof_percentile=None
    ):

        if ground_percentile is None:
            ground_percentile = self.paramaters["GroundPercentile"]
        if roof_percentile is None:
            roof_percentile = self.paramaters["RoofPercentile"]

        self._builder_cm = _pybuilder.ComputeBuildingHeights(
            self._builder_cm, dtm, ground_percentile, roof_percentile
        )

    def to_JSON(self, outfile):
        _pybuilder.WriteCityModelJSON(self._builder_cm, outfile)

    def from_JSON(self, infile):
        self._builder_pc = _pybuilder.ReadCityModelJSON(infile)

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
