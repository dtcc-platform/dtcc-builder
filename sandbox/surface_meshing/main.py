import dtcc_builder
from dtcc_model import Mesh, Building, Surface, MultiSurface
from dtcc_builder.meshing import mesh_multisurface, mesh_surface
import dtcc_viewer

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
ms = MultiSurface(surfaces=[surface1, surface2])

mesh = mesh_multisurface(ms, max_mesh_edge_size=5)
mesh.view()
