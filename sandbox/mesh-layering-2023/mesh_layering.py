from build import _mesh_layering
from dtcc import *
from dtcc_builder.model import builder_volume_mesh_to_volume_mesh
from dtcc_builder.logging import info
import numpy as np

def crop_mesh(mesh, x_max: float = 0.05, y_max: float = 0.05):
    
    num_vertices = len(mesh.vertices)
    num_faces = len(mesh.faces)
    info(f"Loaded Mesh with: \n- Num of vertices:\t{num_vertices}\n- Num of faces: \t{num_faces}")
    elevation = mesh.vertices[:,2]
    domain_min_z = np.min(elevation)
    domain_max_z = np.max(elevation)
    domain_min_x = np.min(mesh.vertices[:,0])
    domain_max_x = np.max(mesh.vertices[:,0])
    domain_min_y = np.min(mesh.vertices[:,1])
    domain_max_y = np.max(mesh.vertices[:,1])
    info("Base Surface Mesh Domain")
    info(f"X Range: [{domain_min_x}, {domain_max_x}]")
    info(f"Y Range: [{domain_min_y}, {domain_max_y}]")
    info(f"Z Range: [{domain_min_z}, {domain_max_z}]")
    
    
    new_max_x = x_max * domain_max_x
    new_max_y=  y_max * domain_max_y
    mesh_tile = Mesh()
    vertices = []
    vertices_indexes = []
    for i,v in enumerate(mesh.vertices):
        if v[0] <= new_max_x and v[1] <= new_max_y:
            print(i,v)
            vertices.append(v)
            vertices_indexes.append(i)
    faces = []
    face_indexes = []
    tile_face_colors = []
    for i,face in enumerate(mesh.faces):
        if all(v in vertices_indexes for v in face):
            faces.append(face)
            face_indexes.append(i)
            # tile_face_colors.append(face_colors[i])
    old_faces = faces.copy()
    for i,idx in enumerate(vertices_indexes):
        for j,f in enumerate(faces):
            faces[j] = [i if vertex == idx else vertex for vertex in f]
    for of, nf in zip(old_faces,faces):
        of = np.array([mesh.vertices[of[0]], mesh.vertices[of[1]],mesh.vertices[of[2]]])
        nf = np.array([vertices[nf[0]], vertices[nf[1]], vertices[nf[2]]])
        if not np.array_equal(of,nf):
            print("Error:", of ,nf)
            
    mesh_tile.vertices = np.array(vertices)
    mesh_tile.faces = np.array(faces)
    return mesh_tile

def main():
     
    mesh = load_mesh("../ansys-test-case-2023/output/surface_mesh_tc_1.pb")
    # mesh = crop_mesh(mesh)
    info(mesh.__str__())
    
    _builder_volume_mesh = _mesh_layering.create_volume_mesh(mesh.vertices,mesh.faces, mesh.normals.flatten())
    info("Layering Complete!")
    volume_mesh = builder_volume_mesh_to_volume_mesh(_builder_volume_mesh)
    info("Converted from C++ to python DTCC mesh")
    # volume_mesh.save("output/test1.vtu")
    info(volume_mesh.__str__())
    
if __name__ == "__main__":
    main()