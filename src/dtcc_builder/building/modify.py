from dtcc_model import Building, GeometryType, MultiSurface, Surface
from dtcc_builder.polygons.polygons import (
    polygon_merger,
    simplify_polygon,
    remove_slivers,
    fix_clearance,
)
from dtcc_builder.register import register_model_method
from shapely.geometry import Polygon
from dtcc_builder.logging import debug, info, warning, error

from typing import List, Tuple
from statistics import mean
import numpy as np


@register_model_method
def get_footprint(building: Building, geom_type: GeometryType = None) -> Surface:
    geom = None
    if geom_type is not None:
        geom = building.flatten_geometry(geom_type)
    else:
        lod_levels = [
            GeometryType.LOD0,
            GeometryType.LOD1,
            GeometryType.LOD2,
            GeometryType.LOD3,
        ]
        for lod in lod_levels:
            geom = building.flatten_geometry(lod)
            if geom is not None:
                break

    if geom is None:
        warning(f"Building {building.id} has no LOD geometry.")
        return None
    height = geom.bounds.zmax
    footprint = geom.to_polygon()
    # print(f"get footprint has {len(footprint.exterior.coords)} vertices")
    s = Surface()
    s.from_polygon(footprint, height)
    return s


def merge_building_footprints(
    buildings: List[Building],
    lod: GeometryType = None,
    max_distance: float = 0.5,
    min_area=10,
) -> List[Building]:
    buildings_geom = [building.flatten_geometry(lod) for building in buildings]
    # print(buildings_geom)
    buildings_geom = [geom for geom in buildings_geom if geom is not None]
    building_heights = [geom.zmax for geom in buildings_geom]
    footprints = [geom.to_polygon(max_distance / 4) for geom in buildings_geom]
    merged_footprint, merged_indices = polygon_merger(
        footprints, max_distance, min_area=min_area
    )
    merged_buildings = []
    for idx, footprint in enumerate(merged_footprint):
        if footprint.geom_type == "MultiPolygon":
            ValueError("Merged footprint is a MultiPolygon")
        footprint = footprint.simplify(1e-2, True)
        if footprint.geom_type == "MultiPolygon":
            ValueError("simplified footprint is a MultiPolygon")
        footprint = remove_slivers(footprint, max_distance / 2)
        if footprint.geom_type == "MultiPolygon":
            ValueError("de-slivered footprint is a MultiPolygon")
        indices = merged_indices[idx]
        height = sum([building_heights[i] * footprints[i].area for i in indices]) / sum(
            [footprints[i].area for i in indices]
        )

        building_surface = Surface()
        building_surface.from_polygon(footprint, height)

        building = Building()
        building.add_geometry(building_surface, GeometryType.LOD0)
        merged_buildings.append(building)
    return merged_buildings
