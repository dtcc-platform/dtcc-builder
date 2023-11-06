import dtcc_io as io
import dtcc_builder

# import dtcc_viewer

# footprints = "tests/data/MinimalCase/PropertyMap.shp"
footprints = "data/helsingborg-residential-2022/footprints.shp"

# las_file = "tests/data/MinimalCase/pointcloud.las"
las_file = "data/helsingborg-residential-2022/pointcloud.las"

city = io.load_city(footprints)
pc = io.load_pointcloud(las_file)

_, bounds = dtcc_builder.builders.calculate_bounds(footprints, las_file)


p = dtcc_builder.parameters.default()
p["min_building_distance"] = 1.0
p["mesh_resolution"] = 1.0
p["ground_smoothing"] = 0
p["auto_domain"] = False
p["x_min"] = 0
p["y_min"] = 0


city = dtcc_builder.builders.build_city(city, pc, bounds, p)
terrain_mesh = dtcc_builder.meshing.terrain_mesh(city, 1.0)
surface_mesh = dtcc_builder.builders.build_city_surface_mesh(city, p, True)
io.save_mesh(surface_mesh, "surface_mesh.obj")
# terrain_mesh.view(pc=pc)
# surface_mesh.view()
