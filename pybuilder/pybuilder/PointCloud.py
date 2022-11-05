import _pybuilder
import os
import numpy


def readLasFiles(las_path, extra_data=True):
    las_path = str(las_path)

    # print(f"loading las from {las_path}")
    if os.path.isdir(las_path):
        pc = _pybuilder.LASReadDirectory(las_path, extra_data)
    elif os.path.isfile(las_path):
        pc = _pybuilder.LASReadFile(las_path, extra_data)
    return pc


def getLasBounds(las_path):
    las_path = str(las_path)
    if not os.path.isdir(las_path):
        print(f"{las_path} is not a directory")
        return None
    bbox = _pybuilder.LASBounds(las_path)
    return bbox


def globalOutlierRemover(point_cloud, outlierMargin):
    point_cloud = _pybuilder.GlobalOutlierRemover(point_cloud, outlierMargin)
    return point_cloud


def VegetationFilter(point_cloud):
    point_cloud = _pybuilder.VegetationFilter(point_cloud)
    return point_cloud


def pointCloud2numpy(point_cloud):
    pts = point_cloud.points
    pts = [[p.x, p.y, p.z] for p in pts]
    return numpy.array(pts)
