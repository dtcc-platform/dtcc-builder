from dtcc_model import Mesh, VolumeMesh
import numpy as np
from math import sqrt
from time import time



def timer(func):
    # This function shows the execution time of
    # the function object passed
    def wrap_func(*args, **kwargs):
        t1 = time()
        result = func(*args, **kwargs)
        t2 = time()
        print(f"Function {func.__name__!r} executed in {(t2-t1):.4f}s")
        return result

    return wrap_func


def distance(point_a, point_b):
    squared_distance = (
        (point_a[0] - point_b[0]) ** 2
        + (point_a[1] - point_b[1]) ** 2
        + (point_a[2] - point_b[2]) ** 2
    )
    return sqrt(squared_distance)


def unit_vector(vector):
    """Returns the unit vector of the vector."""
    return vector / distance(vector, (0, 0, 0))
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """Returns the angle in radians between vectors 'v1' and 'v2'::

    >>> angle_between((1, 0, 0), (0, 1, 0))
    1.5707963267948966
    >>> angle_between((1, 0, 0), (1, 0, 0))
    0.0
    >>> angle_between((1, 0, 0), (-1, 0, 0))
    3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def get_vertex_adj_faces(mesh: Mesh):
    vertex_neighbours = [[] for _ in mesh.vertices]

    for i, face in enumerate(mesh.faces):
        vertex_neighbours[face[0]].append(i)
        vertex_neighbours[face[1]].append(i)
        vertex_neighbours[face[2]].append(i)
    return vertex_neighbours

def get_vertex_adj_cells(volume_mesh: VolumeMesh):
    vertex_neighbours = [[] for _ in volume_mesh.vertices]

    for c, cell in enumerate(volume_mesh.cells.astype(int)):
        vertex_neighbours[cell[0]].append(c)
        vertex_neighbours[cell[1]].append(c)
        vertex_neighbours[cell[2]].append(c)
        vertex_neighbours[cell[3]].append(c)
    return vertex_neighbours




@timer
def get_face_adj_faces(mesh: Mesh):
 
    # Create a dictionary to store neighbor cells for each cell
    neighbor_faces = {i: [] for i in range(mesh.num_faces)}

    # Create a dictionary to store faces of each cell
    face_edges = {i: [] for i in range(mesh.num_faces)}
    shared_edges = {i: [] for i in range(mesh.num_faces)}
    # Populate cell_faces with the faces of each cell
    # for i, face in enumerate(mesh.faces):
    #     # Generate all possible faces by removing one vertex at a time
    #     for j in range(3):
    #         edge = np.delete(face, j)
    #         face_edges[i].append(tuple(sorted(edge)))

    # Iterate through each pair of cells to find neighbors
    for i in range(mesh.num_faces):
        for j in range(i + 1, mesh.num_faces):
            if len(neighbor_faces[i]) == 3: 
                break
            # Check if the two cells share exactly three vertices
            shared_vertices = np.intersect1d(mesh.faces[i], mesh.faces[j])
            if len(shared_vertices) == 2:
                # # Check if the shared vertices form a common face for both cells
                # shared_face_i = cell_faces[i][cell_faces[i].index(tuple(sorted(shared_vertices)))]
                # shared_face_j = cell_faces[j][cell_faces[j].index(tuple(sorted(shared_vertices)))]
                # if shared_face_i == shared_face_j:
                #     neighbor_cells[i].append(j)
                #     neighbor_cells[j].append(i)
                neighbor_faces[i].append(j)
                shared_edges[i].append(shared_vertices)
                
                neighbor_faces[j].append(i)
                shared_edges[j].append(shared_vertices)
            
    return neighbor_faces,shared_edges

@timer
def get_cell_adj_cells(volume_mesh):
 
    # Create a dictionary to store neighbor cells for each cell
    neighbor_cells = {i: [] for i in range(volume_mesh.num_cells)}

    # Create a dictionary to store faces of each cell
    cell_faces = {i: [] for i in range(volume_mesh.num_cells)}
    shared_faces = {i: [] for i in range(volume_mesh.num_cells)}
    # Populate cell_faces with the faces of each cell
    # for i, cell in enumerate(volume_mesh.cells):
    #     # Generate all possible faces by removing one vertex at a time
    #     for j in range(4):
    #         face = np.delete(cell, j)
    #         cell_faces[i].append(tuple(sorted(face)))

    # Iterate through each pair of cells to find neighbors
    for i in range(volume_mesh.num_cells):
        for j in range(i + 1, volume_mesh.num_cells):
            if len(neighbor_cells[i]) == 4: 
                break
            # Check if the two cells share exactly three vertices
            shared_vertices = np.intersect1d(volume_mesh.cells[i], volume_mesh.cells[j])
            if len(shared_vertices) == 3:
                # # Check if the shared vertices form a common face for both cells
                # shared_face_i = cell_faces[i][cell_faces[i].index(tuple(sorted(shared_vertices)))]
                # shared_face_j = cell_faces[j][cell_faces[j].index(tuple(sorted(shared_vertices)))]
                # if shared_face_i == shared_face_j:
                #     neighbor_cells[i].append(j)
                #     neighbor_cells[j].append(i)
                neighbor_cells[i].append(j)
                shared_faces[i].append(shared_vertices)
                
                neighbor_cells[j].append(i)
                shared_faces[j].append(shared_vertices)
    return neighbor_cells,shared_faces