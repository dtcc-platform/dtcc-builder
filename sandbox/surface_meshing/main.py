import dtcc_builder
from dtcc_model import Mesh, Building, Surface, MultiSurface, GeometryType
from dtcc_builder.meshing import mesh_multisurface, mesh_surface
import dtcc_viewer
import dtcc_io
from time import time

surface1 = Surface(
    vertices=[
        [0, 0, 5],
        [0, 10, 5],
        [10, 10, 8],
        [10, 0, 8],
    ]
)

surface2 = Surface(
    vertices=[[0, 0, 0], [10, 0, 0], [10, 10, 0], [2, 10, 0], [2, 12, 0], [0, 12, 0]]
)
# mesh = mesh_surface(surface1, max_mesh_edge_size=5)

start_time = time()
city = dtcc_io.load_cityjson(
    "../../../dtcc-io/sandbox/cityjson/data/DA13_3D_Buildings_Merged.city.json"
)

print(f"Load city took {time() - start_time} seconds")
buildings = city.buildings
print(f"Number of buildings: {len(buildings)}")
start_time = time()
building_multi_surfaces = [b.flatten_geometry(GeometryType.LOD2) for b in buildings]
print(f"Extracting multi surfaces took {time() - start_time} seconds")
print(f"Number of multi surfaces: {len(building_multi_surfaces)}")

start_time = time()
naive_mesh = dtcc_builder.meshing.mesh_multisurfaces(
    building_multi_surfaces, max_mesh_edge_size=-20
)
print(f"Naive meshing took {time() - start_time} seconds")


start_time = time()
merged_mesh = dtcc_builder.meshing.merge_meshes(naive_mesh)
print(f"Merging took {time() - start_time} seconds")

print(f"Number of vertices: {len(merged_mesh.vertices)}")
print(f"Number of faces: {len(merged_mesh.faces)}")
merged_mesh.view()
# start_time = time()
# fancy_mesh = dtcc_builder.meshing.mesh_multisurfaces(
#     building_multi_surfaces, max_mesh_edge_size=5
# )
# print(f"Fancy meshing took {time() - start_time} seconds")
# mesh.view()
