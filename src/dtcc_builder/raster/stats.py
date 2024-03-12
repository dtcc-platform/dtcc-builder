import rasterstats
from dtcc_model import Raster
from shapely.geometry import Polygon
from typing import Union, List
from affine import Affine
from dtcc_builder.register import register_model_method


@register_model_method
def stats(raster: Raster, polygons: Union[Polygon, List[Polygon]], stats=["mean"]):
    """
    Compute statistics for a raster within a polygon.

    Args:
        polygons: A Polygon or a list of Polygons.
        stats: A list of statistics to compute. Supported statistics are:
        ['count', 'min', 'max', 'mean', 'median', 'majority', 'minority', 'unique',  'sum', 'std', 'var', 'percentile_X' (where X is a number between 0 and 100)]

    """
    if isinstance(polygons, Polygon):
        polygons = [polygons]

    stats_str = " ".join(stats)
    rstats = rasterstats.zonal_stats(
        polygons, raster.data, affine=raster.georef, stats=stats_str
    )
    if len(stats) == 1:
        rstats = [s[stats_str] for s in rstats]
    if len(polygons) == 1:
        rstats = rstats[0]
    return rstats
