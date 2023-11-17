import dtcc_builder
import dtcc_io

from pathlib import Path
import numpy as np

in_dir = Path(__file__).parent / ".." / ".." / "data" / "helsingborg-residential-2022"
out_dir = Path(__file__).parent / "output"
out_dir.mkdir(exist_ok=True)

city = dtcc_io.load_city(in_dir / "footprints.shp")
pc = dtcc_io.load_pointcloud(in_dir / "pointcloud.las")

city = dtcc_builder.build_city(city, pc, bounds=None)

offset = (-city.bounds.xmin, -city.bounds.ymin, 0)


ground_mesh = dtcc_builder.meshing.terrain_mesh(city, 5.0, 30.0, 3, False)
ground_mesh.translate(*offset)
ground_mesh.save(out_dir / "ground_mesh.stl")

buildings = dtcc_builder.meshing.building_meshes(
    city, 5.0, 30.0, False, True, False, True
)
buildings = buildings.translate(*offset)
buildings.save(out_dir / "buildings.stl")

merged_city = city.merge_buildings(2, 15, height_merge_strategy="area_weighted")


surface_mesh = dtcc_builder.meshing.city_surface_mesh(
    merged_city, 5, 30, 3, merge_meshes=True, merge_buildings=False
)
surface_mesh = surface_mesh.translate(*offset)
surface_mesh.save(out_dir / "surface_mesh.stl")
