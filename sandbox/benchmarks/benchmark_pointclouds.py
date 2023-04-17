import dtcc_io as io
from dtcc_model.pointcloud import Pointcloud
from time import time
import numpy as np

from dtcc_builder import _pybuilder, pointcloud_processing

start_time = time()
pc = io.load_pointcloud(
    "/home/cit-wastberg/dtcc-platform/dtcc-builder/data/HelsingborgResidential2022/PointCloud.las"
)
print(f"Loaded {len(pc)} points in {time() - start_time} seconds")


def py_outliers(pc, margin):
    z_pts = pc.points[:, 2]
    z_mean = np.mean(z_pts)
    z_std = np.std(z_pts)
    return np.where(np.abs(z_pts - z_mean) > margin * z_std)[0]


start_time = time()
b_pc = _pybuilder.createBuilderPointCloud(
    pc.points, pc.classification, pc.return_number, pc.number_of_returns
)
print(f"Created builder pointcloud in {time() - start_time} seconds")

start_time = time()
outliers = py_outliers(pc, 1)
pc.remove_points(outliers)
print(f"Found {len(outliers)} outliers in {time() - start_time} seconds")

start_time = time()
b_pc = _pybuilder.createBuilderPointCloud(
    pc.points, pc.classification, pc.return_number, pc.number_of_returns
)
print(f"Created builder pointcloud in {time() - start_time} seconds")

start_time = time()
b_pc = _pybuilder.GlobalOutlierRemover(b_pc, 1)
print(f"Found outliers in {time() - start_time} seconds")

start_time = time()
veg_removed = pointcloud_processing.remove_vegetation(pc)
print(f"Removed vegetation in {time() - start_time} seconds")

start_time = time()
b_vef_removed = _pybuilder.VegetationFilter(b_pc)
print(f"Removed vegetation pybind in {time() - start_time} seconds")
