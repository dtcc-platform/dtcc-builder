import dtcc_builder
import dtcc_io

from pathlib import Path
import numpy as np

from triMeshQuality import TriMeshQuality


in_dir = Path(__file__).parent / ".." / ".." / "data" / "helsingborg-residential-2022"
out_dir = Path(__file__).parent / "output"
out_dir.mkdir(exist_ok=True)

city = dtcc_io.load_city(in_dir / "footprints.shp")
pc = dtcc_io.load_pointcloud(in_dir / "pointcloud.las")

city = dtcc_builder.build_city(city, pc, bounds=None)

offset = (-city.bounds.xmin, -city.bounds.ymin, 0)


# ground_mesh = dtcc_builder.meshing.terrain_mesh(city, 5.0, 30.0, 3, False)
# ground_mesh.translate(*offset)
# ground_mesh.save(out_dir / "ground_mesh.stl")

# buildings = dtcc_builder.meshing.building_meshes(
#     city, 5.0, 30.0, False, True, False, True
# )
# buildings = buildings.translate(*offset)
# buildings.save(out_dir / "buildings.stl")

merged_city = city.merge_buildings(2, 15, height_merge_strategy="area_weighted")
merged_city = merged_city.simplify_buildings(0.1)

results = {}

for tc in [0.1, 0.5, 1, 2, 3, 5]:
    fixed_city = merged_city.fix_building_clearance(tc, 15)
    fixed_city.save(out_dir / f"fixed_city_tc_{tc}.gpkg")
    # surface_mesh = dtcc_builder.meshing.city_surface_mesh(
    #     fixed_city, 5, 30, 3, merge_meshes=True, merge_buildings=False
    # )
    surface_mesh = dtcc_builder.meshing.terrain_mesh(fixed_city, 3, 30, 3, True)
    surface_mesh = surface_mesh.translate(*offset)
    aspect_ratio = TriMeshQuality.aspect_ratio(surface_mesh)
    skew = TriMeshQuality.skewness(surface_mesh)
    area = TriMeshQuality.area(surface_mesh)

    results[tc] = {
        "aspect_ratio": aspect_ratio,
        "skew": skew,
        "area": area,
        "num_faces": surface_mesh.faces.shape[0],
    }
    surface_mesh.save(out_dir / f"terrain_mesh_tc_{tc}.stl")

for k, v in results.items():
    print(f"tc: {k}")
    print(
        f"aspect_ratio: {v['aspect_ratio'].min()}, {np.percentile(v['aspect_ratio'], 1)}, {np.percentile(v['aspect_ratio'], 99)}, {v['aspect_ratio'].max()}"
    )
    print(
        f"skew: {v['skew'].min()}, {np.percentile(v['skew'], 1)} {np.percentile(v['skew'], 99)} {v['skew'].max()}"
    )
    print(
        f"area: {v['area'].min()}, {np.percentile(v['area'], 1)}, {np.percentile(v['area'], 99)}, {v['area'].max()}"
    )
# surface_mesh = dtcc_builder.meshing.city_surface_mesh(
#     merged_city, 5, 30, 3, merge_meshes=True, merge_buildings=False
# )
# surface_mesh = surface_mesh.translate(*offset)


# surface_mesh.save(out_dir / "surface_mesh.stl")
