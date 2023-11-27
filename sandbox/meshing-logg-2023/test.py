from dtcc import *
from pathlib import Path
from tetraMeshQuality import check_volume_mesh
from triMeshQuality import check_surface_mesh
import os
import shapely

# Set data paths
data_directory = Path("data/helsingborg-residential-2022")
buildings_path = data_directory / "footprints.shp"
pointcloud_path = data_directory

# Set parameters
p = parameters.default()
p["output_directory"] = "output"
p["verbose"] = True
p["debug"] = True
p["auto_domain"] = False
x0 = 102000.0
y0 = 6213000.0
p["x0"] = x0
p["y0"] = y0
p["x_min"] = 200.0
p["y_min"] = 350.0
p["x_max"] = 400.0
p["y_max"] = 500.0
p["min_building_detail"] = 1.5
p["max_mesh_size"] = 10.0
p["min_mesh_angle"] = 30.0


# Calculate bounds
origin, bounds = calculate_bounds(buildings_path, pointcloud_path, p)
print(bounds)

# Load data from file
city = load_city(buildings_path, bounds=bounds)
pointcloud = load_pointcloud(pointcloud_path, bounds=bounds)

# Build city model
city = build_city(city, pointcloud, bounds, p)
city.view()
city = city.merge_buildings(p["min_building_detail"], p["min_building_area"])

city = city.simplify_buildings(p["min_building_detail"] / 4)
# city.view()
city = city.remove_small_buildings(p["min_building_area"])
city.view()
city = city.fix_building_clearance(p["min_building_detail"], p["min_building_angle"])
city.view()


exit()
# Build city mesh and volume mesh (tetrahedral mesh)
vmesh, bmesh = build_volume_mesh(city, p)

# Save to file
if not os.path.exists("output"):
    os.makedirs("output")
vmesh.save("output/vmesh.pb")
bmesh.save("output/bmesh.pb")
vmesh.save("output/vmesh.vtu")
bmesh.save("output/bmesh.vtu")

exit()

# Compute mesh quality
check_volume_mesh(vmesh)
check_surface_mesh(bmesh)
