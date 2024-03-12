from dtcc_model import Mesh, Surface, MultiSurface, Mesh
from register import register_model_method
import _dtcc_builder
from dtcc_builder.model import (
    create_builder_multisurface,
    create_builder_surface,
    builder_mesh_to_mesh,
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
    if triangle_size is None or triangle_size <= 0:
        triangle_size = -1
    builder_mesh = _dtcc_builder.mesh_multisurface(builder_ms, triangle_size, weld=weld)
    mesh = builder_mesh_to_mesh(builder_mesh)
    return mesh


def mesh_surface(s: Surface, triangle_size=None, weld=False) -> Mesh:
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
    builder_mesh = _dtcc_builder.mesh_surface(builder_surface, triangle_size, weld=weld)
    mesh = builder_mesh_to_mesh(builder_mesh)
    return mesh
