from dtcc import *
from pathlib import Path
from tetraMeshQuality import check_volume_mesh
from triMeshQuality import check_surface_mesh
import os


# Set data paths
data_directory = Path("data/helsingborg-residential-2022")
buildings_path = data_directory / "footprints.shp"
pointcloud_path = data_directory

# Set parameters
p = parameters.default()
p["auto_domain"] = False
x0 = 102000.0
y0 = 6213000.0
p["x0"] = x0
p["y0"] = y0
p["x_min"] = 0.0
p["y_min"] = 0.0
p["x_max"] = 50.0
p["y_max"] = 50.0
p["min_vertex_distance"] = 5.0
p["mesh_resolution"] = 5.0

# Calculate bounds
origin, bounds = calculate_bounds(buildings_path, pointcloud_path, p)
print(bounds)

# Load data from file
city = load_city(buildings_path, bounds=bounds)
pointcloud = load_pointcloud(pointcloud_path, bounds=bounds)

# Build city model
city = build_city(city, pointcloud, bounds, p)

# Build city mesh and volume mesh (tetrahedral mesh)
vmesh, bmesh = build_volume_mesh(city, p)

# Save to file
if not os.path.exists("output"):
    os.makedirs("output")
vmesh.save("output/vmesh.pb")
bmesh.save("output/bmesh.pb")
vmesh.save("output/vmesh.vtu")
bmesh.save("output/bmesh.vtu")

# Compute mesh quality
check_volume_mesh(vmesh)
check_surface_mesh(bmesh)
