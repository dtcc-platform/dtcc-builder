from dtcc_model import Mesh, Surface, MultiSurface
from dtcc_builder.register import register_model_method

from dtcc_builder import _dtcc_builder

from dtcc_builder.model import (
    create_builder_multisurface,
    create_builder_surface,
    builder_mesh_to_mesh,
    mesh_to_builder_mesh,
)


def mesh_multisurface(ms: MultiSurface, triangle_size=None, weld=False) -> Mesh:
    """
    Mesh a `MultiSurface` object into a `Mesh` object.

    Args:
        triangle_size (float): The maximum size of the triangles in the mesh (default None, no max size).
        weld (bool): Whether to weld the vertices of the mesh (default False).

    Returns:
        Mesh: A `Mesh` object representing the meshed `MultiSurface`.
    """

    builder_ms = create_builder_multisurface(ms)
    min_mesh_angle = 25
    if triangle_size is None or triangle_size <= 0:
        triangle_size = -1
    builder_mesh = _dtcc_builder.mesh_multisurface(
        builder_ms, triangle_size, min_mesh_angle, weld
    )
    mesh = builder_mesh_to_mesh(builder_mesh)
    return mesh


def mesh_surface(s: Surface, triangle_size=None) -> Mesh:
    """
    Mesh a `Surface` object into a `Mesh` object.

    Args:
        triangle_size (float): The maximum size of the triangles in the mesh (default None, no max size).
        weld (bool): Whether to weld the vertices of the mesh (default False).

    Returns:
        Mesh: A `Mesh` object representing the meshed `Surface`.
    """

    builder_surface = create_builder_surface(s)
    if triangle_size is None or triangle_size <= 0:
        triangle_size = -1
    builder_mesh = _dtcc_builder.mesh_surface(builder_surface, triangle_size, 25)
    mesh = builder_mesh_to_mesh(builder_mesh)
    return mesh


def mesh_multisurfaces(
    multisurfaces: [MultiSurface], max_mesh_edge_size=-1, min_mesh_angle=25, weld=False
) -> [Mesh]:
    # start_time = time()
    # print(f"flatten multisurfaces took {time() - start_time} seconds")
    # start_time = time()
    builder_multisurfaces = [create_builder_multisurface(ms) for ms in multisurfaces]
    # print(f"create builder multisurfaces took {time() - start_time} seconds")
    # start_time = time()
    meshes = _dtcc_builder.mesh_multisurfaces(
        builder_multisurfaces, max_mesh_edge_size, min_mesh_angle, weld
    )
    # print(f"mesh multisurfaces took {time() - start_time} seconds")
    # start_time = time()
    meshes = [builder_mesh_to_mesh(mesh) for mesh in meshes]
    # print(f"convert builder mesh to mesh took {time() - start_time} seconds")
    return meshes


def merge_meshes(meshes: [Mesh], weld=False) -> Mesh:
    builder_meshes = [mesh_to_builder_mesh(mesh) for mesh in meshes]
    merged_mesh = _dtcc_builder.merge_meshes(builder_meshes, weld)
    mesh = builder_mesh_to_mesh(merged_mesh)
    return mesh


@register_model_method
def merge(mesh: Mesh, other: Mesh, weld=False) -> Mesh:
    builder_mesh = mesh_to_builder_mesh(mesh)
    builder_other = mesh_to_builder_mesh(other)
    merged_mesh = _dtcc_builder.merge_meshes([builder_mesh, builder_other], weld)
    mesh = builder_mesh_to_mesh(merged_mesh)
    return mesh
