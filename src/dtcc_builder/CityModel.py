from google.protobuf.json_format import MessageToJson

from dtcc_builder import _pybuilder

import numpy

import dtcc_io as io
from dtcc_builder import PointCloud, Parameters, ElevationModel
import dtcc_model as model

# from dtccpybuilder.PointCloud import PointCloud
# from dtccpybuilder.ElevationModel import ElevationModel
from typing import List, Tuple

# from dtccpybuilder.Parameters import load_parameters


class CityModel:
    def __init__(
        self, footprints, parmeters = None):
        if parameters is None:
            parameters = Parameters.load_parameters()
        self.parameters = parameters
        self._builder_cm = None
        self.origin = (0, 0)
        self.cleaned = False
        self.simplified = False
        self.extracted_points = False
        self.calculated_heights = False
        
        self.load_protobuf(footprints)

    def get_buildings(self):
        if self._builder_cm is None:
            return None
        else:
            return self._builder_cm.buildings

    def set_buildings(self, b):
        return None

    buildings = property(get_buildings)

    # def generate_citymodel(self, shp_footprint_file):
    #     """load shp file of building footprints"""
    #     if self.bounds is None:
    #         self.bounds = io.citymodel.building_bounds(
    #             shp_footprint_file, self.parameters["DomainMargin"]
    #         )
    #     self._builder_cm = _pybuilder.GenerateCityModel(
    #         str(shp_footprint_file),
    #         self.bounds,
    #         float(self.parameters["MinBuildingDistance"]),
    #         float(self.parameters["MinBuildingSize"]),
    #     )

    def set_origin(self, origin: Tuple[float, float]):
        """set the origin for the city model. Everything will be offset so that origin is at (0,0)"""
        self._builder_cm = _pybuilder.SetCityModelOrigin(self._builder_cm, origin)
        self.origin = origin

    def clean_citymodel(self, min_vert_distance=None):
        """fix any errors in the building polygons and simplify complex polygons"""
        if min_vert_distance is None:
            min_vert_distance = self.parameters["MinVertexDistance"]
        self._builder_cm = _pybuilder.CleanCityModel(
            self._builder_cm, min_vert_distance
        )
        self.cleaned = True

    def simplify_citymodel(self, bounds):
        self._builder_cm = _pybuilder.SimplifyCityModel(
            self._builder_cm,
            bounds,
            self.parameters["MinBuildingDistance"],
            self.parameters["MinVertexDistance"],
        )
        self.simplified = True

    def extract_building_points(
        self,
        point_cloud: PointCloud,
        ground_margins=None,
        ground_outlier_margin=None,
    ):
        """Extract roof and ground points from the point cloud for each building footprint"""

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
        """use RANSAC to filter extreme outliers from roof points. Only use this on very noisy and unclassified data. On clean data it can make things worse"""

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
        """use statistical outlier method to filter roof points for each building"""

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
        """compute the height of each building based on the roof and ground points.
        building_points_statistical_outlier_remover must have been called before calling this method"""
        if ground_percentile is None:
            ground_percentile = self.parameters["GroundPercentile"]
        if roof_percentile is None:
            roof_percentile = self.parameters["RoofPercentile"]

        self._builder_cm = _pybuilder.ComputeBuildingHeights(
            self._builder_cm, dtm._grid_field, ground_percentile, roof_percentile
        )
        self.calculated_heights = True

    def get_building_roof_points(self, building_idx):
        roof_pts = self.buildings[building_idx].roof_points
        roof_pts = [[p.x, p.y, p.z] for p in roof_pts]
        return numpy.array(roof_pts)

    def to_JSON(self, outfile):
        """serialize CItyModel to a JSON file"""
        pbcm = model.CityModel()
        with open(outfile, "w") as f:
            f.write(MessageToJson(pbcm.FromString(self.to_protobuf())))

    # def from_JSON(self, infile):
    #     """Load CityModel from JSON file"""
    #     self._builder_cm = _pybuilder.ReadCityModelJSON(str(infile))

    def load_protobuf(self, protobuf: model.CityModel):
        """load CityModel from a CityModel protobuf string"""
        self._builder_cm = _pybuilder.loadCityModelProtobuf(protobuf.SerializeToString())
        self.origin = (self._builder_cm.origin.x, self._builder_cm.origin.y)
        
    
    def to_protobuf(self) -> str:
        """convert CityModel to protobuf string"""
        pb_string = _pybuilder.convertCityModelToProtobuf(self._builder_cm)
        return pb_string


def get_building_footprint(building: _pybuilder.Building):
    footprint = building.footprint
