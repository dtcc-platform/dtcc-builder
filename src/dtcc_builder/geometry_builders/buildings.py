from dtcc_model import Building, GeometryType
from dtcc_model import Surface, MultiSurface, PointCloud, Raster
from dtcc_builder.model import create_builder_polygon, create_builder_pointcloud
from dtcc_builder import _dtcc_builder
from dtcc_builder.logging import debug, info, warning, error
from shapely.geometry import Polygon
import numpy as np

from typing import List, Tuple

from .surface import extrude_surface


def extrude_building(
    building: Building, default_ground_height=0, always_use_default=False
) -> MultiSurface:
    """
    Extrudes the LOD0 representation of a building from its height to the ground leve.

    Parameters
    ----------
    `building` : dtcc_model.Building
        The building to extrude.
    `default_ground_height` : float, optional
        If building does not have a ground_height property, the default ground
        level to use, by default 0.
    `always_use_default_ground` : bool, optional
        Whether to always use the default ground height or use ground_height attribute, by default False.

    Returns
    -------
    `MultiSurface`
        The extruded building.
    """
    if always_use_default:
        ground_height = default_ground_height
    else:
        ground_height = building.attributes.get("ground_height", default_ground_height)

    geometry = building.lod0
    if geometry is None:
        error(f"Building {building.id} has no LOD0 geometry.")
        return None
    if isinstance(geometry, Surface):
        geometry = MultiSurface(surfaces=[geometry])
    if not isinstance(geometry, MultiSurface):
        error(f"Building {building.id} LOD0 geometry is not a (Multi)Surface.")
        return None

    extrusion = MultiSurface()

    for surface in geometry.surfaces:
        extrusion = extrusion.merge(extrude_surface(surface, ground_height))
    return extrusion


def compute_building_heights(
    buildings: List[Building],
    terrain: Raster,
    min_building_height=2.5,
    roof_percentile=0.9,
    overwrite=False,
) -> List[Building]:
    info("Computing building heights...")
    for building in buildings:
        footprint = building.lod0
        if footprint is None:
            warning(f"Building {building.id} has no LOD0 geometry.")
            continue
        centroid = footprint.centroid
        ground_height = terrain.get_value(centroid[0], centroid[1])
        building.attributes["ground_height"] = ground_height
        if overwrite or footprint.zmax == 0:
            roof_points = building.point_cloud
            if roof_points is None or len(roof_points) == 0:
                warning(f"Building {building.id} has no roof points. using min height")
                footprint.set_z(ground_height + min_building_height)
            else:
                z_values = roof_points.points[:, 2]
                roof_top = np.percentile(z_values, roof_percentile * 100)
                height = roof_top - ground_height
                if height < min_building_height:
                    height = min_building_height
                footprint.set_z(ground_height + height)
    return buildings


def build_lod1_buildings(
    buildings: [Building],
    default_ground_height=0,
    always_use_default_ground=False,
    rebuild=True,
) -> List[Building]:
    """
    Build the LOD1 representation of the given buildings.

    Parameters
    ----------
    `buildings` : [dtcc_model.Building]
        The buildings to build the LOD1 representation of.
    `default_ground_height` : float, optional
        If building does not have a ground_height property, the default ground
        level to use, by default 0.
    `always_use_default_ground` : bool, optional
        Whether to always use the default ground height or use groun_height attribute, by default False.
    `rebuild` : bool, optional
        Whether to rebuild the LOD1 representation if it already exists, by default True.
    Returns
    -------
    [dtcc_model.Building]
        The buildings with the LOD1 representation built.
    """
    info(f"Building LOD1 representations of {len(buildings )} buildings...")
    for building in buildings:
        if building.lod1 is not None and not rebuild:
            continue
        if building.lod0 is None:
            warning(f"Building {building.id} has no LOD0 geometry.")
            continue
        geometry = extrude_building(
            building, default_ground_height, always_use_default_ground
        )
        if geometry is not None:
            building.add_geometry(geometry, GeometryType.LOD1)
        else:
            warning(f"Building {building.id} LOD1 geometry could not be built.")
    return buildings


def extract_roof_points(
    buildings: [Building],
    pointcloud: PointCloud,
    statistical_outlier_remover=True,
    roof_outlier_neighbors=5,
    roof_outlier_margin=1.5,
    ransac_outlier_remover=False,
    ransac_outlier_margin=3.0,
    ransac_iterations=250,
) -> List[Building]:
    footprint_polygons = [b.get_footprint() for b in buildings]

    builder_polygon = [
        create_builder_polygon(p.to_polygon())
        for p in footprint_polygons
        if p is not None
    ]
    if len(pointcloud.points) == len(pointcloud.classification):
        ground_mask = np.logical_or(
            pointcloud.classification == 2, pointcloud.classification == 9
        )
        not_ground_mask = ~ground_mask
        points = pointcloud.points[not_ground_mask]
    else:
        points = pointcloud.points
    roof_points = _dtcc_builder.extract_building_points(
        builder_polygon,
        points,
        statistical_outlier_remover,
        roof_outlier_neighbors,
        roof_outlier_margin,
    )

    idx = 0
    # some buildings may not have a footprint, and thus not have roof points
    for fp in footprint_polygons:
        if fp is not None:
            pc = PointCloud(points=roof_points[idx])
            pc.calculate_bounds()
            buildings[idx].add_geometry(pc, GeometryType.POINT_CLOUD)
            idx += 1
    return buildings
