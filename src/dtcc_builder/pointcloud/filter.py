import numpy as np
from typing import List

from dtcc_model import PointCloud
from dtcc_builder.logging import info, warning, error
from dtcc_builder.register import register_model_method


@register_model_method
def remove_global_outliers(pc: PointCloud, margin: float):
    """
    Remove outliers from a `PointCloud` whose Z-values are more than `margin`
    standard deviations from the mean.

    Args:
        margin (float): The margin in standard deviations to consider a point an outlier.

    Returns:
        PointCloud: A new `PointCloud` object with the outliers removed.
    """

    z_pts = pc.points[:, 2]
    z_mean = np.mean(z_pts)
    z_std = np.std(z_pts)
    outliers = np.where(np.abs(z_pts - z_mean) > margin * z_std)[0]
    new_pc = pc.copy()
    new_pc.remove_points(outliers)
    return new_pc


@register_model_method
def classification_filter(pc: PointCloud, classes: List[int], keep: bool = False):
    """
    Filter a `PointCloud` object based on its classification.

    Args:
        classes (List[int]): The classification values to keep or remove.
        keep (bool): Whether to keep the points with the specified classification values (default False, remove them).

    Returns:
        PointCloud: A new `PointCloud` object with the specified points removed.
    """
    if len(pc.points) != len(pc.classification):
        warning("Pointcloud not classified, returning original pointcloud.")
        return pc
    mask = np.isin(pc.classification, classes)
    if keep:
        mask = np.logical_not(mask)
    pc.remove_points(mask)
    return pc


@register_model_method
def remove_vegetation(pc: PointCloud) -> PointCloud:
    """
    Return a pioint cloud with vegetation removed.

    Args:
        pc (PointCloud): The `PointCloud` object to remove vegetation from.

    Returns:
        PointCloud: A new `PointCloud` object with the vegetation removed.
    """
    new_pc = pc.copy()
    veg_indices = _find_vegetation(pc)
    new_pc.remove_points(veg_indices)
    return new_pc


def _find_vegetation(pc: PointCloud, filter_on_return_number=True):
    """Find the indices of points that belong to vegetation in a point cloud.

    Args:
        pc: A `PointCloud` object representing the point cloud to filter.
        filter_on_return_number: A boolean indicating whether to filter on return number (default True).

    Returns:
        A 1D NumPy array of indices of points that belong to vegetation.
    """

    has_classification = len(pc.classification) == len(pc.points)
    has_return_number = len(pc.return_number) == len(pc.points)
    if not has_classification and not has_return_number:
        warning(
            "Classification and return number are not set for all points. Ignoring vegetation filter."
        )
        return np.array([])
    if not has_classification:
        warning("Classification is not set for all points. Ignoring")

    if filter_on_return_number and not has_return_number:
        filter_on_return_number = False
        warning("Return number is not set for all points. Ignoring")

    classes_with_vegetation = set([3, 4, 5])
    used_classes = pc.used_classifications()
    veg_classes = classes_with_vegetation.intersection(used_classes)
    if len(veg_classes) == 0:
        has_classification = False
    else:
        veg_classes = np.array(list(veg_classes))
        filter_on_return_number = False

    vegetation_indices = np.array([])
    if has_classification:
        vegetation_indices = np.where(np.isin(pc.classification, veg_classes))[0]

    elif filter_on_return_number:
        is_veg = pc.return_number != pc.num_returns

        # only reclassify points that are not already classified
        if len(pc.classification) == len(pc.points):
            is_veg = np.logical_and(is_veg, pc.classification == 1)
        vegetation_indices = np.where(is_veg)[0]

    return vegetation_indices
