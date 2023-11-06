from dtcc import *
from pathlib import Path
from volume_mesh_quality import check_volume_mesh
from surface_mesh_quality import check_surface_mesh

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
vmesh.save("vmesh.pb")
bmesh.save("bmesh.pb")
vmesh.save("vmesh.vtu")
bmesh.save("bmesh.vtu")

# Compute mesh quality
check_volume_mesh(vmesh)
check_surface_mesh(bmesh)
