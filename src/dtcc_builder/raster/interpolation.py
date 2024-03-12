import scipy.ndimage
import numpy as np
import rasterio
from dtcc_builder.register import register_model_method
from dtcc_model import Raster
from logging import info, warning, error


@register_model_method
def fill_holes(raster: Raster):
    """
    Fill nodata holes in a `Raster` object using the nearest neighbour.

    Returns:
        Raster: A new `Raster` object with the holes filled.
    """
    data = raster.data
    nodata = raster.nodata
    mask = data == nodata
    if np.any(mask):
        info(f"filling {mask.sum()} holes in raster")
        ind = scipy.ndimage.distance_transform_edt(
            mask, return_distances=False, return_indices=True
        )
        data = data[tuple(ind)]
    raster.data = data
    return raster


@register_model_method
def resample(raster: Raster, cell_size=None, scale=None, method="bilinear"):
    """
    Resample a `Raster` object to a new cell size or scale.

    Args:
        cell_size (float): The new cell size in meters (default None).
        scale (float): The scaling factor for the cell size (default None).
        method (str): The resampling method to use, one of "bilinear", "nearest", or "cubic" (default "bilinear").

    Returns:
        Raster: A new `Raster` object with the resampled data.
    """
    sample_methods = {
        "bilinear": 1,
        "nearest": 0,
        "cubic": 3,
    }
    if cell_size is None and scale is None:
        raise ValueError("Either cell_size or scale must be specified")
    if not method in sample_methods:
        raise ValueError(
            f"Invalid resampling method, use one of {list(sample_methods.keys())}"
        )
    _raster = raster.copy()
    if cell_size is not None:
        scale = cell_size / _raster.cell_size[0]
    if scale == 1:
        return _raster
    _raster.data = scipy.ndimage.zoom(
        _raster.data,
        scale,
        order=sample_methods[method],
        mode="nearest",
        grid_mode=True,
    )
    _raster.georef *= _raster.georef.scale(scale, scale)

    return _raster
