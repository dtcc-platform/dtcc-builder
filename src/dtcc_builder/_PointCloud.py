from google.protobuf.json_format import MessageToJson

from dtcc_builder import _dtcc_builder
import dtcc_io as io
import dtcc_model as model
import os

from pathlib import Path

import numpy
from typing import List, Tuple


class PointCloud:
    """
    Class for storing point cloud object.

    This class encapsulates functionality for working with point cloud data,
    providing methods for loading, processing, and converting point cloud data
    in various formats. It allows loading point clouds from LAS files, converting
    to/from protobuf format, performing various point cloud processing tasks,
    and exporting point clouds to different formats.

    Attributes
    ----------
    `origin` : Tuple[float, float]
        Origin coordinates of the point cloud.
    `bounds` : tuple
        Bounding box of the point cloud.
    `used_classifications` : dict
        Dictionary of used point classifications.

    Parameters
    ----------
    `pointcloud_path` : str or None, optional
        Path to a point cloud file (LAS format), by default None.
    `pb_sting` : bytes or None, optional
        Protobuf string representing the point cloud, by default None.
    `bounds` : tuple, optional
        Bounding box coordinates, by default ().

    Raises
    ------
    Exception
        If point cloud cannot be loaded from the given path.
    """

    def __init__(self, pointcloud_path=None, pb_string=None, bounds=()):
        self._builder_pc = None
        self.origin = (0, 0)
        self.bounds = bounds
        self.used_classifications = {}
        if pointcloud_path is not None:
            self.load_from_path(pointcloud_path)
        if pb_string is not None:
            self.from_protobuf(pb_string)

    def __len__(self):
        """
        Return the number of points in the point cloud.

        Returns
        -------
        `int`
            Number of points in the point cloud.
        """
        if self._builder_pc is None:
            return 0
        else:
            return len(self._builder_pc)

    def _points_as_numpy(self):
        """
        Convert internal points to a NumPy array.

        Returns
        -------
        `numpy.ndarray`
            NumPy array containing point coordinates.
        """
        pts = self._builder_pc.points
        pts = [[p.x, p.y, p.z] for p in pts]
        return numpy.array(pts)

    points = property(_points_as_numpy)
    """
    NumPy array containing point coordinates.

    Returns
    -------
    `numpy.ndarray`
        NumPy array containing point coordinates.
    """

    def load_from_path(self, las_path, extra_data=True):
        """
        Load point cloud from a LAS file and set internal representation.

        Parameters
        ----------
        `las_path` : str
            Path to the LAS file.
        `extra_data` : bool, optional
            Flag indicating whether to load extra data, by default True.

        Raises
        ------
        Exception
            If point cloud cannot be loaded from the given path.
        """
        pb_string = io.load_pointcloud(
            las_path,
            points_classification_only=not extra_data,
            bounds=self.bounds,
            return_serialized=True,
        )
        if pb_string is not None:
            self.from_protobuf(pb_string)
        else:
            raise Exception(f"Could not load point cloud from path {las_path}")

    def set_origin(self, origin: Tuple[float, float]):
        """
        Set the origin coordinates of the point cloud.

        Parameters
        ----------
        `origin` : tuple of float
            Origin coordinates as (x, y).
        """
        self.origin = origin
        self._builder_pc = _dtcc_builder.SetPointCloudOrigin(self._builder_pc, origin)

    def remove_global_outliers(self, outlier_margin):
        """
        Remove global outliers from the point cloud.

        Parameters
        ----------
        `outlier_margin` : float
            Margin for outlier removal.
        """
        self._builder_pc = _dtcc_builder.GlobalOutlierRemover(
            self._builder_pc, outlier_margin
        )

    def vegetation_filter(self):
        """
        Apply a vegetation filter to the point cloud.
        """
        self._builder_pc = _dtcc_builder.VegetationFilter(self._builder_pc)

    def get_bounds(self):
        """
        Get the bounding box of the point cloud.

        Returns
        -------
        `tuple`
            Bounding box coordinates (x_min, y_min, x_max, y_max).
        """
        pts = pointcloud_to_numpy(self._builder_pc)
        x_min, y_min, z_min = pts.min(axis=0)
        x_max, y_max, z_max = pts.max(axis=0)
        return (x_min, y_min, x_max, y_max)

    def from_protobuf(self, protobuf_string: bytes):
        """
        Load point cloud from a protobuf string.

        Parameters
        ----------
        `protobuf_string` : bytes
            Protobuf string representing the point cloud.
        """
        self._builder_pc = _dtcc_builder.loadPointCloudProtobuf(protobuf_string)
        self.used_classification = set(self._builder_pc.classifications)

    def to_protobuf(self) -> str:
        """
        Convert the point cloud to a protobuf string.

        Returns
        -------
        `str`
            Protobuf string representing the point cloud.
        """
        return _dtcc_builder.convertPointCloudToProtobuf(self._builder_pc)

    def to_json(self, output_path: Path):
        """
        Convert the point cloud to JSON format and save it to a file.

        Parameters
        ----------
        `output_path` : Path
            Path to save the JSON file.
        """
        pbpc = model.PointCloud()
        with open(output_path, "w") as f:
            f.write(MessageToJson(pbpc.FromString(self.to_protobuf())))
