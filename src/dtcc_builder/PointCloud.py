from google.protobuf.json_format import MessageToJson

from dtcc_builder import _pybuilder
import dtcc_io as io
import os

from pathlib import Path

import numpy
from typing import List, Tuple


class PointCloud:
    """Class for storing point cloud object"""

    def __init__(self, pointcloud_path=None, pb_sting = None, bounds=()):
        """ """
        self._builder_pc = None
        self.origin = (0, 0)
        self.bounds = bounds
        self.used_classifications = {}
        if pointcloud_path is not None:
            self.load_from_path(pointcloud_path)
        if pb_sting is not None:
            self.from_protobuf(pb_sting)

    def __len__(self):
        if self._builder_pc is None:
            return 0
        else:
            return len(self._builder_pc)

    def _points_as_numpy(self):
        pts = self._builder_pc.points
        pts = [[p.x, p.y, p.z] for p in pts]
        return numpy.array(pts)

    points = property(_points_as_numpy)

    def load_from_path(self, las_path, extra_data=True):
        pb_string = io.pointcloud.read(las_path, points_classification_only = not extra_data, bounds = self.bounds, return_serialized=True)
        if pb_string is not None:
            self.from_protobuf(pb_string)
        else:
            raise Exception(f"Could not load point cloud from path {las_path}")

    def set_origin(self, origin: Tuple[float, float]):
        self.origin = origin
        self._builder_pc = _pybuilder.SetPointCloudOrigin(self._builder_pc, origin)

    def remove_global_outliers(self, outlier_margin):
        self._builder_pc = _pybuilder.GlobalOutlierRemover(
            self._builder_pc, outlier_margin
        )

    def vegetation_filter(self):
        self._builder_pc = _pybuilder.VegetationFilter(self._builder_pc)

    def get_bounds(self):
        pts = pointcloud_to_numpy(self._builder_pc)
        x_min, y_min, z_min = pts.min(axis=0)
        x_max, y_max, z_max = pts.max(axis=0)
        return (x_min, y_min, x_max, y_max)

    def from_protobuf(self, protobuf_string: bytes):
        self._builder_pc = _pybuilder.loadPointCloudProtobuf(protobuf_string)
        self.used_classification = set(self._builder_pc.classifications)

    def to_protobuf(self) -> str:
        return _pybuilder.convertPointCloudToProtobuf(self._builder_pc)

    def to_json(self, output_path: Path):
        pbpc = io.dtcc_model.protobuf.dtcc_pb2.PointCloud()
        with open(output_path, "w") as f:
            f.write(MessageToJson(pbpc.FromString(self.to_protobuf())))
        
