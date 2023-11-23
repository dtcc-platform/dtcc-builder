import dtcc_builder
import dtcc_io

from pathlib import Path
import numpy as np

from triMeshQuality import TriMeshQuality

in_dir = Path(__file__).parent / ".." / ".." / "data" / "helsingborg-residential-2022"
out_dir = Path(__file__).parent / "output"
out_dir.mkdir(exist_ok=True)

city = dtcc_io.load_city(in_dir / "footprints.shp")

merged_city = city.merge_buildings(2, 15, height_merge_strategy="area_weighted")
merged_city = merged_city.simplify_buildings(0.1)
target_clearance = 2.0
footprints = [b.footprint for b in city.buildings]
    for idx, f in enumerate(footprints):

