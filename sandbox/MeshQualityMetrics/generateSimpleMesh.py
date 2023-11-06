import numpy as np
import dtcc_io
from dtcc_model import Mesh, VolumeMesh



def generate_unit_box(nx, ny) -> Mesh:
    '''The function generates a mesh representing a unit box with a given number of subdivisions in the x
    and y directions.
    
    Parameters
    ----------
    nx
        The parameter `nx` represents the number of subdivisions along the x-axis of the unit box.
    ny
        The parameter "ny" represents the number of divisions along the y-axis of the unit box. It
    determines the number of vertices and faces that will be generated in the y-direction.
    
    Returns
    -------
        a Mesh object.
    
    '''

    hx = 1 / nx
    hy = 1 / ny

    V = []
    K = []

    for iy in range(ny + 1):
        for ix in range(nx + 1):
            x = ix*hx
            y = iy*hy
            V.append((x, y, 0.0))

    for iy in range(ny):
        for ix in range(nx):
            v0 = iy*(nx + 1) + ix
            v1 = v0 + 1
            v2 = v0 + (nx + 1)
            v3 = v1 + (nx + 1)
            K.append((v0, v1, v3))
            K.append((v0, v3, v2))

    mesh = Mesh()
    mesh.vertices = np.array(V)
    mesh.faces = np.array(K)

    return mesh

def generate_unit_cube(nx: int, ny: int, nz: int) -> VolumeMesh:
    '''The function `generate_unit_cube` generates a volume mesh of a unit cube with specified dimensions.
    
    Parameters
    ----------
    nx : int
        The parameter `nx` represents the number of subdivisions along the x-axis of the unit cube.
    ny : int
        The parameter `ny` represents the number of divisions along the y-axis in the unit cube. It
    determines the number of vertices and cells that will be generated in the y-direction.
    nz : int
        The parameter `nz` represents the number of divisions along the z-axis in the unit cube. It
    determines the number of vertices and cells that will be generated in the z-direction.
    
    Returns
    -------
        a VolumeMesh object.
    
    '''

    hx = 1 / nx
    hy = 1 / ny
    hz = 1 / nz 
    
    V = []
    Cells = []

    for iz in range(nz+1):
        for iy in range(ny + 1):
            for ix in range(nx + 1):
                x = ix*hx
                y = iy*hy
                z = iz*hz
                V.append((x, y, z))
    
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                v0 = iz*(ny + 1)*(nx + 1) +iy*(nx + 1) + ix
                v1 = v0 + 1
                v2 = v0 + (nx + 1)
                v3 = v1 + (nx + 1)
                v4 = v0 + (ny + 1)* (nx + 1)
                v5 = v1 + (ny + 1)* (nx + 1)
                v6 = v2 + (ny + 1)* (nx + 1)
                v7 = v3 + (ny + 1)* (nx + 1)
                
                # Create tetrahedron cells
                Cells.append((v0, v1, v3, v7))
                Cells.append((v0, v1, v5, v7))
                Cells.append((v0, v4, v5, v7))
                Cells.append((v0, v2, v3, v7))
                Cells.append((v0, v4, v6, v7))
                Cells.append((v0, v2, v6, v7))

    mesh = VolumeMesh()
    mesh.vertices = np.array(V)
    mesh.cells = np.array(Cells)

    return mesh

def generate_cube(size, resolution) -> VolumeMesh:

    
    nx = resolution[0]
    ny = resolution[1]
    nz = resolution[2] 
    
    hx = size[0] / resolution[0]
    hy = size[1] / resolution[1]
    hz = size[2] / resolution[2] 
    
    V = []
    Cells = []

    for iz in range(nz+1):
        for iy in range(ny + 1):
            for ix in range(nx + 1):
                x = ix*hx
                y = iy*hy
                z = iz*hz
                V.append((x, y, z))
    
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                v0 = iz*(ny + 1)*(nx + 1) +iy*(nx + 1) + ix
                v1 = v0 + 1
                v2 = v0 + (nx + 1)
                v3 = v1 + (nx + 1)
                v4 = v0 + (ny + 1)* (nx + 1)
                v5 = v1 + (ny + 1)* (nx + 1)
                v6 = v2 + (ny + 1)* (nx + 1)
                v7 = v3 + (ny + 1)* (nx + 1)
                
                # Create tetrahedron cells
                Cells.append((v0, v1, v3, v7))
                Cells.append((v0, v1, v5, v7))
                Cells.append((v0, v4, v5, v7))
                Cells.append((v0, v2, v3, v7))
                Cells.append((v0, v4, v6, v7))
                Cells.append((v0, v2, v6, v7))

    mesh = VolumeMesh()
    mesh.vertices = np.array(V)
    mesh.cells = np.array(Cells)

    return mesh





