from dtccpybuilder import _pybuilder
import os
import numpy
from typing import List, Tuple


def calc_las_bounds(las_path):
    las_path = str(las_path)
    if not os.path.isdir(las_path):
        print(f"{las_path} is not a directory")
        return None
    bbox = _pybuilder.LASBounds(las_path)
    return bbox


class PointCloud:
    def __init__(self, las_path=None, bounds=()):
        self._builder_pc = None
        self.origin = (0, 0)
        self.bounds = bounds
        self.used_classifications = {}
        if las_path is not None:
            self.read_las_files(las_path)

    def __len__(self):
        if self._builder_pc is None:
            return 0
        else:
            return len(self._builder_pc)

    def points_as_numpy(self):
        pts = self._builder_pc.points
        pts = [[p.x, p.y, p.z] for p in pts]
        return numpy.array(pts)

    points = property(points_as_numpy)

    def read_las_files(self, las_path, extra_data=True):
        las_path = str(las_path)

        # print(f"loading las from {las_path}")
        if os.path.isdir(las_path):
            pc = _pybuilder.LASReadDirectory(las_path, self.bounds, extra_data)
        elif os.path.isfile(las_path):
            pc = _pybuilder.LASReadFile(las_path, self.bounds, extra_data)
        self._builder_pc = pc
        self.used_classification = set(self._builder_pc.classifications)

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

    def init_with_protobuf(self, protobuf_string: str):
        self._builder_pc = _pybuilder.loadPointCloudProtobuf(protobuf_string)

    def to_protobuf(self) -> str:
        return _pybuilder.convertPointCloudToProtobuf(self._builder_pc)
