# Test mesh coloring based on dyadic mesh sizes
# Anders Logg 2023-11-27

from dtcc import *

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

# Test mesh from ansys-test-case-2023 sandbox (run first)
mesh = load_mesh("../ansys-test-case-2023/output/surface_mesh_tc_1.pb")
num_vertices = len(mesh.vertices)
num_faces = len(mesh.faces)
#mesh.view()

def ideal_layer_height(area):
    # Compute ideal layer height for regular tetrahedron
    c = 2.0**1.5*3**(-0.75)
    h = c*np.sqrt(area)
    return h

def check_neighbors(colors, ff):
    # Check coloring of neighboring faces
    max_diff = 0
    num_big_diffs = 0
    for i, faces in enumerate(ff):
        diffs = [abs(colors[i] - colors[j]) for j in faces]
        big_diffs = [d for d in diffs if d > 1]
        max_diff = max(max_diff, max(diffs))
        if len(big_diffs) > 0:
            #print(i, diffs)
            num_big_diffs += 1
    print("Max diff:", max_diff)
    print(f"Num big diffs: {num_big_diffs} / {num_faces} ({100*num_big_diffs/num_faces:.2f}%)")

# Compute mesh sizes
areas = np.zeros(num_faces)
for i, f in enumerate(mesh.faces):
    v = [mesh.vertices[j][:2] for j in f]
    a = v[1] - v[0]
    b = v[2] - v[0]
    A = 0.5*abs(np.cross(a, b))
    areas[i] = A

# Compute ideal layer heights for smallest and largest mesh sizes
_min_height = ideal_layer_height(np.min(areas))
_max_height = ideal_layer_height(np.max(areas))
print("Ideal min layer height:", _min_height)
print("Ideal max layer height:", _max_height)

# Compute dyadic mesh sizes to match min/max as close as possible
rho = _max_height / _min_height
mid = np.sqrt(_min_height*_max_height)
num_layers = int(np.log2(rho) + 0.5)
min_height = mid / 2**(num_layers/2)
layer_heights = np.array([min_height*2**i for i in range(num_layers+1)])
print("Number of layers:", num_layers)
print("Adjusted min layer height:", layer_heights[0])
print("Adjusted max layer height:", layer_heights[-1])
print("Layer heights:", layer_heights)

# Assign layer heights to mesh (closest by quotient)
colors = np.zeros(num_faces, dtype=np.int32)
for i in range(len(areas)):
    h = ideal_layer_height(areas[i])
    d = np.abs(np.log(h / layer_heights))
    color = np.argmin(d)
    colors[i] = color

# Build mapping from vertices to faces
vf = [set() for i in range(num_vertices)]
for i, f in enumerate(mesh.faces):
    for j in f:
        vf[j].add(i)

# Build mapping from faces to faces
ff = [set() for i in range(num_faces)]
for i, f in enumerate(mesh.faces):
    for j in f:
        ff[i] |= vf[j]

# Check coloring of neighboring faces
check_neighbors(colors, ff)

# Plot mesh colors
_mesh = tri.Triangulation(mesh.vertices[:, 0], mesh.vertices[:, 1], mesh.faces)
plt.tripcolor(_mesh, facecolors=colors, edgecolors="k")
#plt.show()
