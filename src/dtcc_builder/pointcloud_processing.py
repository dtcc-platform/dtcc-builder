# Copyright (C) 2022 Dag Wästberg
# Licensed under the MIT License
#
# Modified by Anders Logg 2023
#
# This module provides functionality for point cloud processing.

import numpy as np
import logging


def remove_gloabl_outliers(pc, margin):
    outliers = find_global_outliers(pc, margin)
    pc.remove_points(outliers)
    return pc


def find_global_outliers(pc, margin):
    z_pts = pc.points[:, 2]
    z_mean = np.mean(z_pts)
    z_std = np.std(z_pts)
    return np.where(np.abs(z_pts - z_mean) > margin * z_std)[0]


def remove_vegetation(pc, filter_on_return_number=True):
    vegetation = find_vegetation(pc, filter_on_return_number)
    pc.remove_points(vegetation)
    return pc


def find_vegetation(pc, filter_on_return_number=True):
    has_classification = len(pc.classification) == len(pc.points)
    has_return_number = len(pc.return_number) == len(pc.points)
    if not has_classification and not has_return_number:
        logging.warning(
            "Classification and return number are not set for all points. Ignoring vegetation filter."
        )
        return np.array([])
    if not has_classification:
        logging.warning("Classification is not set for all points. Ignoring")

    if filter_on_return_number and not has_return_number:
        logging.warning("Return number is not set for all points. Ignoring")

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
        vegetation_indices = np.where(pc.return_number != pc.number_of_returns)[0]

    return vegetation_indices
