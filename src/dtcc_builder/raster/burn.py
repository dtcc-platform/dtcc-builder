import rasterio
from rasterio.features import rasterize

from dtcc_model import Raster
from dtcc_builder.register import register_model_method
from shapely.geometry import Polygon


@register_model_method
def burn_polygons(
    raster: Raster, poly_values: tuple, all_touched: bool = True
) -> Raster:
    """
    Burn polygons into a raster.

    Args:
        raster (Raster): The input raster.
        poly_values (tuple): A tuple of (polygon, value) pairs to burn into the raster.
        all_touched (bool, optional): Whether to rasterize all pixels touched by the polygon, or just those whose center is within the polygon. Defaults to True.

    Returns:
        Raster: A new raster with the polygons burned in.
    """
    out_raster = raster.copy()
    arr = out_raster.data

    transform = rasterio.transform.from_bounds(
        raster.bounds.west,
        raster.bounds.south,
        raster.bounds.east,
        raster.bounds.north,
        raster.width,
        raster.height,
    )
    for poly, value in poly_values:
        mask = rasterize(
            [poly],
            out_shape=arr.shape,
            transform=transform,
            fill=0,
            all_touched=all_touched,
        )
        arr[mask == 1] = value

    out_raster.data = arr
    return out_raster
