# Copyright (C) 2023 Anders Logg
# Licensed under the MIT License
#
# This demo illustrates how to build a city model from raw data,
# essentially equivalent to running the dtcc-build command-line
# utility, but with more control over the process.

from dtcc import *
from pathlib import Path

# Set data paths
data_directory = Path("data/miljovis2023")
buildings_path = data_directory / "by_04.shp"
pointcloud_path = data_directory

# Set parameters
p = parameters.default()
p["auto_domain"] = True

# Calculate bounds
origin, bounds = calculate_bounds(buildings_path, pointcloud_path, p)

# Load data from file
city = load_city(buildings_path, bounds=bounds)
pointcloud = load_pointcloud(pointcloud_path, bounds=bounds)

# Build city model
city = build_city(city, pointcloud, bounds, p)

# Build ground mesh and building mesh (surface meshes)
ground_mesh, building_mesh = build_mesh(city, p)

# Build city mesh and volume mesh (tetrahedral mesh)
volume_mesh, volume_mesh_boundary = build_volume_mesh(city, p)

# Save data to file
city.save(data_directory / "city.pb")
ground_mesh.save(data_directory / "ground_mesh.pb")
building_mesh.save(data_directory / "building_mesh.pb")
volume_mesh.save(data_directory / "volume_mesh.pb")
volume_mesh_boundary.save(data_directory / "volume_mesh_boundary.pb")

# View data
city.view()
pointcloud.view()
volume_mesh_boundary.view()
