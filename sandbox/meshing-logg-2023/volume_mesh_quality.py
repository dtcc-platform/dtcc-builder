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

def find_neighbour_cells(volume_mesh: VolumeMesh):
    # Initialize a dictionary to store neighbors for each cell
    print("neighbor_cells!!!")
    neighbor_cells = [[] for _ in range(volume_mesh.num_cells)]

    adj_cells = get_adjacent_cells(volume_mesh)


    for v in range(volume_mesh.num_cells):
        print(v,"/",volume_mesh.num_cells)
        for i in range(len(adj_cells[v])):
            for j in range(i+1,len(adj_cells[v])):
                common_vertices = np.intersect1d(volume_mesh.cells[adj_cells[v][i]],volume_mesh.cells[adj_cells[v][j]])
                if len(common_vertices) == 3:
                    neighbor_cells[i].append(j)
                    neighbor_cells[j].append(i)
    return neighbor_cells

    # # Iterate through all pairs of cells to find neighbors
    # for i in range(volume_mesh.num_cells):
    #     for j in range(i + 1, volume_mesh.num_cells):
    #         common_vertices = cell_vertices[i].intersection(cell_vertices[j])

    #         # If two cells share 3 vertices, they are neighbors
    #         if len(common_vertices) == 3:
    #             neighbor_cells[i].append(j)
    #             neighbor_cells[j].append(i)
    # return neighbor_cells
@timer
def surface_aspect_ratio(mesh: Mesh):
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
def surface_orthogonal_quality(mesh: Mesh)->np.ndarray:
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
def surface_skewness(mesh: Mesh):
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
def surface_smoothness(mesh: Mesh):
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
    mean, max, ar = surface_aspect_ratio(mesh)
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

    skew = surface_skewness(mesh)
    print("MAX\t", np.nanmax(skew))
    print("MIN\t", np.nanmin(skew))
    print("MEAN\t", np.nanmean(skew))
    print("MEDIAN\t", np.nanmedian(skew))

    sm = surface_smoothness(mesh)
    print("MAX\t", np.nanmax(sm))
    print("MIN\t", np.nanmin(sm))
    print("MEAN\t", np.nanmean(sm))
    print("MEDIAN\t", np.nanmedian(sm))


#Edge length check provides a lazy aspect ratio indicator for tetrahedron cells
@timer
def volume_aspect_ratio_1(volume_mesh: VolumeMesh):
    aspect_ratio =np.zeros((volume_mesh.num_cells))
    aspect_ratio.fill(np.NaN)

    for c,cell in enumerate(volume_mesh.cells.astype(int)):
        v0 = volume_mesh.vertices[cell[0]]
        v1 = volume_mesh.vertices[cell[1]]
        v2 = volume_mesh.vertices[cell[2]]
        v3 = volume_mesh.vertices[cell[3]]

        e01 = distance(v1, v0)
        e02 = distance(v2, v0)
        e03 = distance(v3, v0)
        e12 = distance(v2, v1)
        e13 = distance(v3, v1)
        e23 = distance(v3, v2)

        min_length = min((e01,e02,e03,e12,e13,e23))
        max_length = max((e01,e02,e03,e12,e13,e23))

        if min_length > 0:
            aspect_ratio[c] = max_length/min_length

    return aspect_ratio

@timer
def volume_aspect_ratio_2(volume_mesh: VolumeMesh):
    aspect_ratio =np.zeros((volume_mesh.num_cells))
    aspect_ratio.fill(np.NaN)

    for c,cell in enumerate(volume_mesh.cells.astype(int)):
        v0 = volume_mesh.vertices[cell[0]]
        v1 = volume_mesh.vertices[cell[1]]
        v2 = volume_mesh.vertices[cell[2]]
        v3 = volume_mesh.vertices[cell[3]]

        e01 = v1 - v0
        e02 = v2 - v0
        e03 = v3 - v0
        e12 = v2 - v1
        e13 = v3 - v1
        e23 = v3 - v2


        d01 = distance(v1, v0)
        d02 = distance(v2, v0)
        d03 = distance(v3, v0)
        d12 = distance(v2, v1)
        d13 = distance(v3, v1)
        d23 = distance(v3, v2)

        #Volume of cell
        V = 1/6 * abs(np.linalg.det([e01, e02,e03]))

        # Area of Faces cells
        a_0 = 0.5 * np.linalg.norm(np.outer(e01, e02))
        a_1 = 0.5 * np.linalg.norm(np.outer(e01, e03))
        a_2 = 0.5 * np.linalg.norm(np.outer(e12, e13))
        a_3 = 0.5 * np.linalg.norm(np.outer(e02, e03))

        r = 3 * V / (a_0+a_1+a_2+a_3)
        # min_length = min((e01,e02,e03,e12,e13,e23))
        max_length = max((d01,d02,d03,d12,d13,d23))

        if r > 0:
            aspect_ratio[c] = max_length/(2*sqrt(6)*r)

    return aspect_ratio

@timer
def volume_aspect_ratio_3(volume_mesh: VolumeMesh):
    aspect_ratio =np.zeros((volume_mesh.num_cells))
    aspect_ratio.fill(np.NaN)

    for c,cell in enumerate(volume_mesh.cells.astype(int)):
        v0 = volume_mesh.vertices[cell[0]]
        v1 = volume_mesh.vertices[cell[1]]
        v2 = volume_mesh.vertices[cell[2]]
        v3 = volume_mesh.vertices[cell[3]]

        e01 = v1 - v0
        e02 = v2 - v0
        e03 = v3 - v0
        e12 = v2 - v1
        e13 = v3 - v1
        e23 = v3 - v2


        d01 = distance(v1, v0)
        d02 = distance(v2, v0)
        d03 = distance(v3, v0)
        d12 = distance(v2, v1)
        d13 = distance(v3, v1)
        d23 = distance(v3, v2)

        #Volume of cell
        V = 1/6 * abs(np.linalg.det([e01, e02,e03]))

        # Area of Faces cells
        a_0 = 0.5 * np.linalg.norm(np.outer(e01, e02))
        a_1 = 0.5 * np.linalg.norm(np.outer(e01, e03))
        a_2 = 0.5 * np.linalg.norm(np.outer(e12, e13))
        a_3 = 0.5 * np.linalg.norm(np.outer(e02, e03))

        # Inradius
        r = 3 * V / (a_0+a_1+a_2+a_3)

        #Circumradius
        p1 = d01*d23
        p2 = d02*d13
        p3 = d03*d12
        cR = sqrt((p1+p2+p3)*(p1-p2+p3)*(p1+p2-p3)*(-p1+p2+p3))/(24*V)

        if r > 0:
            aspect_ratio[c] = cR/r

    return aspect_ratio



# Opposite edges of regular tetrahedron should be perpendicular of each other.
# So the inner product of two opposite edges could be used as a quality metric
# for the cell's orthogonality
@timer
def volume_orthogonality(volume_mesh: VolumeMesh):
    orthogonality =np.zeros((volume_mesh.num_cells))
    orthogonality.fill(np.NaN)

    for c,cell in enumerate(volume_mesh.cells):
        v0 = volume_mesh.vertices[cell[0]]
        v1 = volume_mesh.vertices[cell[1]]
        v2 = volume_mesh.vertices[cell[2]]
        v3 = volume_mesh.vertices[cell[3]]

        e01 = v1 - v0
        e02 = v2 - v0
        e03 = v3 - v0
        e12 = v2 - v1
        e13 = v3 - v1
        e23 = v3 - v2

        n0 = np.outer(e01,e02)
        n1 = np.outer(e01,e03)
        n2 = np.outer(e02,e03)
        n3 = np.outer(e12,e13)

        orthogonality[c] = sqrt(np.dot(n0,n1)**2 +
                                np.dot(n0,n2)**2 +
                                np.dot(n0,n3)**2 +
                                np.dot(n1,n3)**2 +
                                np.dot(n1,n2)**2 +
                                np.dot(n2,n3)**2 )/6
    return orthogonality

@timer
def volume_orthogonal_quality(volume_mesh: VolumeMesh):
    orthogonality =np.zeros((volume_mesh.num_cells))
    orthogonality.fill(np.NaN)
    # orthogonality = []
    for c,cell in enumerate(volume_mesh.cells.astype(int)):
        v0 = volume_mesh.vertices[cell[0]]
        v1 = volume_mesh.vertices[cell[1]]
        v2 = volume_mesh.vertices[cell[2]]
        v3 = volume_mesh.vertices[cell[3]]

        e01 = v1 - v0
        e02 = v2 - v0
        e03 = v3 - v0
        e12 = v2 - v1
        e13 = v3 - v1
        e23 = v3 - v2

        n0 = unit_vector(np.cross(e01,e02))
        n1 = unit_vector(np.cross(e01,e03))
        n2 = unit_vector(np.cross(e02,e03))
        n3 = unit_vector(np.cross(e12,e13))


        orthogonality[c] = sqrt(np.dot(n0,n1)**2 +
                                np.dot(n0,n2)**2 +
                                np.dot(n0,n3)**2 +
                                np.dot(n1,n3)**2 +
                                np.dot(n1,n2)**2 +
                                np.dot(n2,n3)**2 )/6

    return orthogonality

# Measuring the deviation between the cell volume and the ideal cell volume
@timer
def volume_skewness(volume_mesh: VolumeMesh):
    skewness =np.zeros((volume_mesh.num_cells))
    skewness.fill(np.NaN)

    for c,cell in enumerate(volume_mesh.cells.astype(int)):
        v0 = volume_mesh.vertices[cell[0]]
        v1 = volume_mesh.vertices[cell[1]]
        v2 = volume_mesh.vertices[cell[2]]
        v3 = volume_mesh.vertices[cell[3]]

        e01 = v1 - v0
        e02 = v2 - v0
        e03 = v3 - v0
        # e12 = v2 - v1
        # e13 = v3 - v1
        # e23 = v3 - v2


        #Volume of cell
        V =  abs(np.linalg.det([e01, e02,e03]))/6

        d01 = distance(v1, v0)
        d02 = distance(v2, v0)
        d03 = distance(v3, v0)
        d12 = distance(v2, v1)
        d13 = distance(v3, v1)
        d23 = distance(v3, v2)

        #Circumradius
        p1 = d01*d23
        p2 = d02*d13
        p3 = d03*d12
        R = sqrt((p1+p2+p3)*(p1-p2+p3)*(p1+p2-p3)*(-p1+p2+p3))/(24*V)

        V_ideal = (8*sqrt(3)/27) * (R**3)

        skewness[c] = 1 - V/V_ideal

    return skewness

@timer
def volume_smoothness(volume_mesh: VolumeMesh):
    # smoothness = np.zeros((volume_mesh.num_vertices))
    # smoothness.fill(np.NaN)
    smoothness = []
    neighbouring_cells = get_adjacent_cells(volume_mesh)

    V = np.zeros((volume_mesh.num_cells))
    for c,cell in enumerate(volume_mesh.cells.astype(int)):

        v0 = volume_mesh.vertices[cell[0]]
        v1 = volume_mesh.vertices[cell[1]]
        v2 = volume_mesh.vertices[cell[2]]
        v3 = volume_mesh.vertices[cell[3]]


        e01 = v1 - v0
        e02 = v2 - v0
        e03 = v3 - v0

        #Volume of cell
        V[c] = 1/6 * abs(np.linalg.det([e01, e02,e03]))

    neighbor_cells = find_neighbour_cells(volume_mesh)

    return

def check_volume_mesh(volume_mesh: VolumeMesh):

    print("Smoothness of Volume Mesh:")
    sm = volume_smoothness(volume_mesh)
    print("MAX\t", np.nanmax(sm))
    print("MIN\t", np.nanmin(sm))
    print("MEAN\t", np.nanmean(sm))
    print("MEDIAN\t", np.nanmedian(sm))


def main():

    # mesh_path = "/home/auth-georspai/dtcc-builder/data/HelsingborgResidential2022/bug_output/volume_mesh_boundary.pb"
    # mesh = dtcc_io.load_mesh(mesh_path)
    # check_surface_mesh(mesh)
    # # ar =surface_orthogonal_quality(mesh)
    # # print("MAX\t", np.nanmax(ar))
    # # print("MIN\t", np.nanmin(ar))
    # # print("MEAN\t", np.nanmean(ar))
    # # print("MEDIAN\t", np.nanmedian(ar))


    volume_mesh_path = "/home/auth-georspai/dtcc-builder/data/HelsingborgResidential2022/bug_output/volume_mesh.pb"
    volume_mesh = dtcc_io.load_volume_mesh(volume_mesh_path)


    # ar = volume_aspect_ratio_3(volume_mesh)
    sm = volume_smoothness(volume_mesh)
    #oq = volume_orthogonal_quality(volume_mesh)
    #sk = volume_skewness(volume_mesh)

    # print("Aspect Ratio of Volume Mesh: ")
    # print("MAX\t", np.nanmax(ar))
    # print("MIN\t", np.nanmin(ar))
    # print("MEAN\t", np.nanmean(ar))
    # print("MEDIAN\t", np.nanmedian(ar))

    # print("Skewness of Volume Mesh:")
    # print("MAX\t", np.nanmax(sk))
    # print("MIN\t", np.nanmin(sk))
    # print("MEAN\t", np.nanmean(sk))
    # print("MEDIAN\t", np.nanmedian(sk))

    # print("Orthogonal Quality of Volume Mesh:")
    # print("MAX\t", np.nanmax(oq))
    # print("MIN\t", np.nanmin(oq))
    # print("MEAN\t", np.nanmean(oq))
    # print("MEDIAN\t", np.nanmedian(oq))

    print("Smoothness of Volume Mesh:")
    print("MAX\t", np.nanmax(sm))
    print("MIN\t", np.nanmin(sm))
    print("MEAN\t", np.nanmean(sm))
    print("MEDIAN\t", np.nanmedian(sm))


if __name__ == "__main__":
    main()
