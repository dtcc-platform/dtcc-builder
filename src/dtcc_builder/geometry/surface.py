from dtcc_model import Surface

from register import register_model_method
from meshing.meshing import mesh_surface


@register_model_method
def mesh(s: Surface, triangle_size=None, weld=False) -> Mesh:
    """
    Mesh a `Surface` object into a `Mesh` object.

    Args:
        triangle_size (float): The maximum size of the triangles in the mesh (default None, no max size).
        weld (bool): Whether to weld the vertices of the mesh (default False).

    Returns:
        Mesh: A `Mesh` object representing the meshed `Surface`.
    """

    return mesh_surface(s, triangle_size, weld)
