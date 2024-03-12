from dtcc_builder.polygons.polygons import (
    polygon_merger,
    simplify_polygon,
    remove_slivers,
    fix_clearance,
)

from dtcc_builder.register import register_model_method
from dtcc_model import City, Bounds
from dtcc_model.object.building import Building
from dtcc_model.object.object import GeometryType
from statistics import mean
import shapely
import dataclasses
from copy import deepcopy
from collections import defaultdict
from shapely.geometry import MultiPolygon, Polygon, JOIN_STYLE, CAP_STYLE
from dtcc_builder.logging import info, warning, error


@register_model_method
def simplify_buildings(city: City, tolerance=0.1) -> City:
    """
    Simplify the footprint of buildings in a `City` object.

    Args:
        city (City): The `City` object to simplify the buildings of.
        tolerance (float): The tolerance for simplification (default 0.1).

    Returns:
        City: A new `City` object with the simplified buildings.
    """
    simplified_city = deepcopy(city)
    simplified_city.buildings = []
    for b in city.buildings:
        b = dataclasses.replace(b)
        b.footprint = simplify_polygon(b.footprint, tolerance)
        simplified_city.buildings.append(b)
    return simplified_city


@register_model_method
def remove_small_buildings(city: City, min_area=10) -> City:
    """
    Remove small buildings from a `City` object.

    Args:
        city (City): The `City` object to remove small buildings from.
        min_area (float): The minimum area in square meters for a building to be kept (default 10).

    Returns:
        City: A new `City` object with the small buildings removed.
    """
    filtered_city = deepcopy(city)
    filtered_city.buildings = []
    for b in city.buildings:
        if b.footprint.area > min_area:
            filtered_city.buildings.append(b)
    return filtered_city


@register_model_method
def merge_buildings(
    city: City,
    max_distance=0.15,
    min_area=10,
    simplify=True,
    properties_merge_strategy="list",
    height_merge_strategy="mean",
) -> City:
    """
    Merge buildings that are close together.

    Parameters
    ----------
    max_distance : float, optional
        The maximum distance in meters between buildings to consider them close enough to merge (default 0.15).
    min_area : float, optional
        The minimum area in square meters for a building to be kept (default 10).
    simplify : bool, optional
        Whether to simplify the merged buildings (default True).
    properties_merge_strategy : str, optional
        The strategy for merging properties. Options are 'list' and 'sample'. 'list' will create a list of all properties for the merged building. 'sample' will pick a property value from a random building (default "list").
    height_merge_strategy : str, optional
        The strategy for merging heights. Options are 'mean', 'area_weighted' and 'max' .
        'mean' will take the mean height of the merged buildings.
        'area_weighted' will take the area weighted mean height of the merged buildings.
        'max' will take the maximum height of the merged buildings (default "mean").

    Returns
    -------
    City
        A new `City` object with the merged buildings.
    """
    merged_city = deepcopy(city)
    footprints = [b.footprint for b in city.buildings]
    merged_polygons, merged_polygons_idx = polygon_merger(
        footprints, max_distance, min_area=min_area
    )

    merged_city.buildings = []
    for idx, merged_polygon in enumerate(merged_polygons):
        merged_polygon = merged_polygon.simplify(1e-3, True)
        merged_polygon = remove_slivers(merged_polygon, max_distance / 2)
        if merged_polygon.area < min_area:
            continue
        b = dataclasses.replace(city.buildings[merged_polygons_idx[idx][0]])
        b.footprint = merged_polygon
        if height_merge_strategy == "mean":
            b.height = mean(
                [city.buildings[i].height for i in merged_polygons_idx[idx]]
            )
        elif height_merge_strategy == "area_weighted":
            b.height = sum(
                [
                    city.buildings[i].height * city.buildings[i].footprint.area
                    for i in merged_polygons_idx[idx]
                ]
            ) / sum(
                [city.buildings[i].footprint.area for i in merged_polygons_idx[idx]]
            )
        elif height_merge_strategy == "max":
            b.height = max([city.buildings[i].height for i in merged_polygons_idx[idx]])
        else:
            error(f"Unknown height merge strategy: {height_merge_strategy}")
        b.ground_level = min(
            [city.buildings[i].ground_level for i in merged_polygons_idx[idx]]
        )
        b.roofpoints = city.buildings[merged_polygons_idx[idx][0]].roofpoints
        for i in merged_polygons_idx[idx][1:]:
            b.roofpoints.merge(city.buildings[i].roofpoints)

        property_dicts = [
            city.buildings[i].properties for i in merged_polygons_idx[idx]
        ]
        if properties_merge_strategy == "list":
            merged_properties = defaultdict(list)
            for p in property_dicts:
                for k, v in p.items():
                    merged_properties[k].append(v)
        elif properties_merge_strategy == "sample":
            merged_properties = {}
            for p in property_dicts:
                for k, v in p.items():
                    if v:
                        merged_properties[k] = v
        b.properties = dict(merged_properties)

        merged_city.buildings.append(b)

    return merged_city


@register_model_method
def fix_building_clearance(
    city: City, target_clearance: float, min_angle: float, accepted_tol_fraction=0.9
) -> City:
    """
    Fix the clearance of the footprints in the building models. After running
    each building should have a minimum_clearance of `tol` meters and a minimum
    angle between each edge of `min_angle` degrees.
    """
    fixed_city = city.copy()
    footprints = [b.footprint for b in city.buildings]
    for idx, f in enumerate(footprints):
        fixed_fp = fix_clearance(f, target_clearance, accepted_tol_fraction)
        fixed_city.buildings[idx].footprint = fixed_fp
    fixed_city = fixed_city.merge_buildings(0, 0, height_merge_strategy="area_weighted")

    return fixed_city


# @register_model_method
# def calculate_bounds(city: City, buffer: float = 0) -> City:
#     """
#     Calculate the bounds of a `City` object.

#     Args:
#         city (City): The `City` object to calculate the bounds of.
#         buffer (float): The buffer to add to the bounds (default 0).

#     Returns:
#         City: A new `City` object with the bounds calculated.
#     """
#     footprints = [b.footprint for b in city.buildings]
#     bounds = MultiPolygon(footprints).bounds
#     city.bounds = Bounds(bounds[0], bounds[1], bounds[2], bounds[3])
#     if buffer != 0:
#         city.bounds.buffer(buffer)
#     return city
