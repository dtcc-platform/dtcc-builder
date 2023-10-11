import dtcc_io
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


def get_adjacent_faces(mesh: Mesh):
    vertex_neighbours = [[] for _ in mesh.vertices]

    for i, face in enumerate(mesh.faces):
        vertex_neighbours[face[0]].append(i)
        vertex_neighbours[face[1]].append(i)
        vertex_neighbours[face[2]].append(i)
    return vertex_neighbours

def get_adjacent_cells(volume_mesh: VolumeMesh):
    vertex_neighbours = [[] for _ in volume_mesh.vertices]

    for c, cell in enumerate(volume_mesh.cells.astype(int)):
        vertex_neighbours[cell[0]].append(c)
        vertex_neighbours[cell[1]].append(c)
        vertex_neighbours[cell[2]].append(c)
        vertex_neighbours[cell[3]].append(c)
    return vertex_neighbours


@timer
def aspect_ratio(mesh: Mesh):
    '''The function calculates the average and maximum aspect ratio of the faces in a given mesh.
    
    Parameters
    ----------
    mesh : Mesh
        The `mesh` parameter is an object of the `Mesh` class. It represents a mesh, which is a collection
    of vertices and faces that define the shape of a 3D object. The `Mesh` class likely has attributes
    such as `vertices` (a list of vertex coordinates) and
    
    Returns
    -------
        three values: avg_spaect_ratio (average aspect ratio), max_aspect_ratio (maximum aspect ratio), and
    aspect_ratio (an array containing the aspect ratio for each face of the mesh).
    
    '''
    max_aspect_ratio: float = 0.0
    mean_aspect_ratio: float = 0.0
    aspect_ratio = np.zeros((mesh.num_faces))
    aspect_ratio.fill(np.NaN)
    faces = 0

    for i, (a, b, c) in enumerate(mesh.faces):
        # Calculate the lengths of the edges
        AB = distance(mesh.vertices[a], mesh.vertices[b])
        AC = distance(mesh.vertices[a], mesh.vertices[c])
        BC = distance(mesh.vertices[b], mesh.vertices[c])

        # Find Minimum and Maximum length edge
        l_min = AB
        l_max = AC

        if l_min > BC:
            l_min = BC
        if l_min > AC:
            l_min = AC

        if l_max < BC:
            l_max = BC
        if l_max < AB:
            l_min = AB

        # Calculate the aspect ratio
        if l_min > 0:
            ar = l_max / l_min
            mean_aspect_ratio += ar
            faces = faces + 1
            aspect_ratio[i] = ar
            if ar > max_aspect_ratio:
                max_aspect_ratio = ar

    mean_aspect_ratio = mean_aspect_ratio / faces

    # Some Mesh faces have edges with zero length. these are not taken in consideration.
    print(faces, mesh.num_faces)
    return mean_aspect_ratio, max_aspect_ratio, aspect_ratio


@timer
def surface_orthogonality(mesh: Mesh)->np.ndarray:
    '''The function calculates the surface orthogonality of a mesh by computing the dot products of the
    vectors formed by the mesh vertices and then dividing them by the distances between the vertices.
    
    Parameters
    ----------
    mesh : Mesh
        The `mesh` parameter is an object of type `Mesh`. It represents a mesh structure, typically used in
    computer graphics and computational geometry. The `Mesh` object contains information about the
    vertices and faces of the mesh.
    
    Returns
    -------
        an array of orthogonality values for each face in the given mesh.
    
    '''
    orthogonality = np.zeros((mesh.num_faces))
    orthogonality.fill(np.NaN)

    for i, (a, b, c) in enumerate(mesh.faces):
        AB = mesh.vertices[b] - mesh.vertices[a]
        AC = mesh.vertices[c] - mesh.vertices[a]
        BC = mesh.vertices[c] - mesh.vertices[b]
        
        # Uncomment to prevent 0/0 Runtime Warning
        # if not (np.any(AB) and np.any(BC) and np.any(AC)): continue

        oi_0 = abs(np.dot(AB, AC) / (distance(AB, (0, 0, 0)) * distance(AC, (0, 0, 0))))
        oi_1 = abs(np.dot(AB, BC) / (distance(AB, (0, 0, 0)) * distance(BC, (0, 0, 0))))
        oi_2 = abs(np.dot(AC, BC) / (distance(AC, (0, 0, 0)) * distance(BC, (0, 0, 0))))
        orthogonality[i] = min((oi_0, oi_1, oi_2))
    
    return orthogonality


@timer
def surface_orthogonality_2(mesh: Mesh):
    orthogonality = np.zeros((mesh.num_faces))
    orthogonality.fill(np.NaN)

    for i, (a, b, c) in enumerate(mesh.faces):
        AB = mesh.vertices[b] - mesh.vertices[a]
        AC = mesh.vertices[c] - mesh.vertices[a]
        BC = mesh.vertices[c] - mesh.vertices[b]

        # Uncomment to prevent 0/0 Runtime Warning
        # if not (np.any(AB) and np.any(BC) and np.any(AC)): continue

        oi_0 = abs(1 - (angle_between(AC, AB) / (np.pi / 2)))
        oi_1 = abs(1 - (angle_between(AB, BC) / (np.pi / 2)))
        oi_2 = abs(1 - (angle_between(AC, BC) / (np.pi / 2)))

        orthogonality[i] = min((oi_0, oi_1, oi_2))

    return orthogonality


@timer
def orthogonal_quality(mesh: Mesh)->np.ndarray:
    '''The function calculates the surface orthogonality of a mesh by computing the dot products of the
    vectors formed by the mesh vertices and then dividing them by the distances between the vertices.
    
    Parameters
    ----------
    mesh : Mesh
        The `mesh` parameter is an object of type `Mesh`. It represents a mesh structure, typically used in
    computer graphics and computational geometry. The `Mesh` object contains information about the
    vertices and faces of the mesh.
    
    Returns
    -------
        an array of orthogonality values for each face in the given mesh.
    
    '''
    orthogonality = np.zeros((mesh.num_faces))
    orthogonality.fill(np.NaN)

    for i, (a, b, c) in enumerate(mesh.faces):
        AB = mesh.vertices[b] - mesh.vertices[a]
        AC = mesh.vertices[c] - mesh.vertices[a]
        BC = mesh.vertices[c] - mesh.vertices[b]

        face_normal = mesh.normals[i]
        centroid =(mesh.vertices[a] + mesh.vertices[b] + mesh.vertices[c] ) / 3
        AB_unit_vector = unit_vector(AB)
        AB_normal = np.cross(face_normal, AB_unit_vector) 
        e_AB = unit_vector((mesh.vertices[a] + mesh.vertices[b]) /2 - centroid)
        
        oi_0 =abs( np.dot(e_AB, AB_normal))
        #oi_0 = np.dot(e_AB, AB_normal)
        
        AC_unit_vector = unit_vector(AC)
        AC_normal = np.cross(face_normal, AC_unit_vector)
        e_AC = unit_vector((mesh.vertices[a] + mesh.vertices[c]) /2 - centroid)
        
        oi_1 =abs( np.dot(e_AC, AC_normal))
        #oi_1 =np.dot(e_AC, AC_normal)
        
        BC_unit_vector = unit_vector(BC)
        BC_normal = np.cross(face_normal, BC_unit_vector)
        e_BC = unit_vector((mesh.vertices[b] + mesh.vertices[c]) /2 - centroid)
        
        oi_2 =abs( np.dot(e_BC, BC_normal))
        #oi_2 =np.dot(e_BC, BC_normal)
        # Uncomment to prevent 0/0 Runtime Warning
        # if not (np.any(AB) and np.any(BC) and np.any(AC)): continue

        # oi_0 = abs(np.dot(AB, AC) / (distance(AB, (0, 0, 0)) * distance(AC, (0, 0, 0))))
        # oi_1 = abs(np.dot(AB, BC) / (distance(AB, (0, 0, 0)) * distance(BC, (0, 0, 0))))
        # oi_2 = abs(np.dot(AC, BC) / (distance(AC, (0, 0, 0)) * distance(BC, (0, 0, 0))))
        orthogonality[i] = min((oi_0, oi_1, oi_2))

    return orthogonality


# Equiangular skew
@timer
def skewness(mesh: Mesh):
    '''The function `surface_skewness` calculates the skewness of each face in a given mesh.
    
    Parameters
    ----------
    mesh : Mesh
        The `mesh` parameter is an object of the `Mesh` class. It represents a mesh, which is a collection
    of vertices and faces that define the shape of a 3D object. The `Mesh` class likely has attributes
    such as `vertices` and `faces`, which store the coordinates
    
    Returns
    -------
        an array of skewness values for each face in the given mesh.
    
    '''
    skewness = np.zeros((mesh.num_faces))
    skewness.fill(np.NaN)

    for i, (a, b, c) in enumerate(mesh.faces):
        AB = mesh.vertices[b] - mesh.vertices[a]
        AC = mesh.vertices[c] - mesh.vertices[a]
        BC = mesh.vertices[c] - mesh.vertices[b]

        # Uncomment to prevent 0/0 Runtime Warning
        # if not (np.any(AB) and np.any(BC) and np.any(AC)): continue
        theta_0 = angle_between(AC, AB)
        theta_1 = angle_between(AB, BC)
        theta_2 = angle_between(AC, BC)

        theta_min = min((theta_0, theta_1, theta_2))
        theta_max = max((theta_0, theta_1, theta_2))

        skewness[i] = max(
            (theta_max - np.pi / 3) / (2 * np.pi / 3),
            (np.pi / 3 - theta_min) / (np.pi / 3),
        )

    return skewness


@timer
def smoothness(mesh: Mesh):
    '''The function calculates the smoothness of each vertex in a mesh by computing the standard deviation
    of the areas of the adjacent faces.
    
    Parameters
    ----------
    mesh : Mesh
        The `mesh` parameter is an object of the `Mesh` class. It represents a 3D mesh, which consists of
    vertices and faces. The `Mesh` class likely has attributes such as `vertices` (a list of vertex
    coordinates) and `faces` (a list of vertex indices
    
    Returns
    -------
        an array of smoothness values for each vertex in the mesh.
    
    '''
    smoothness = np.zeros((mesh.num_vertices))
    smoothness.fill(np.NaN)

    neighbouring_faces = get_adjacent_faces(mesh)

    for v, vertex in enumerate(mesh.vertices):
        area = []
        for f in neighbouring_faces[v]:
            AB = mesh.vertices[mesh.faces[f][1]] - mesh.vertices[mesh.faces[f][0]]
            AC = mesh.vertices[mesh.faces[f][2]] - mesh.vertices[mesh.faces[f][0]]
            area.append(0.5 * np.linalg.norm(np.outer(AB, AC)))
        smoothness[v] = np.std(area)

    return smoothness


def check_surface_mesh(mesh: Mesh):
    mean, max, ar = aspect_ratio(mesh)
    print("MAX\t", np.nanmax(ar), max)
    print("MIN\t", np.nanmin(ar))
    print("MEAN\t", np.nanmean(ar), mean)
    print("MEDIAN\t", np.nanmedian(ar))

    c = 0
    for a in ar:
        if a < 5:
            c += 1

    print(c)

    ortho = surface_orthogonality(mesh)
    print("MAX\t", np.nanmax(ortho))
    print("MIN\t", np.nanmin(ortho))
    print("MEAN\t", np.nanmean(ortho))
    print("MEDIAN\t", np.nanmedian(ortho))

    ortho = surface_orthogonality_2(mesh)
    print("MAX\t", np.nanmax(ortho))
    print("MIN\t", np.nanmin(ortho))
    print("MEAN\t", np.nanmean(ortho))
    print("MEDIAN\t", np.nanmedian(ortho))

    skew = skewness(mesh)
    print("MAX\t", np.nanmax(skew))
    print("MIN\t", np.nanmin(skew))
    print("MEAN\t", np.nanmean(skew))
    print("MEDIAN\t", np.nanmedian(skew))

    sm = smoothness(mesh)
    print("MAX\t", np.nanmax(sm))
    print("MIN\t", np.nanmin(sm))
    print("MEAN\t", np.nanmean(sm))
    print("MEDIAN\t", np.nanmedian(sm))




def main():
    path_to_mesh = "/home/auth-georspai/scratch/mariya_meshes/Output/mesh_4_0/volume_mesh_boundary.pb"
    mesh = dtcc_io.load_mesh(path_to_mesh)
    
    skew = skewness(mesh)
    print("MAX\t", np.nanmax(skew))
    print("MIN\t", np.nanmin(skew))
    print("MEAN\t", np.nanmean(skew))
    print("MEDIAN\t", np.nanmedian(skew))


if __name__ == "__main__":
    main()