import dtcc_io
from dtcc_model import Mesh, VolumeMesh

import numpy as np

from qualityMetricsUtils import * 
from triMeshQuality import *




class TetraElementQuality():
    C_tetra = 124.70765802

    @staticmethod
    def volume(v0: np.ndarray,v1: np.ndarray,v2: np.ndarray,v3: np.ndarray)->float:
        return abs(np.linalg.det([v1 - v0, v2 - v0, v3 - v0]))/6
    
    @staticmethod   
    def element_quality(v0: np.ndarray,v1: np.ndarray,v2: np.ndarray,v3: np.ndarray)->float:   
        e01 = v1 - v0
        e02 = v2 - v0
        e03 = v3 - v0 
                
        l01 = distance(v1, v0)
        l02 = distance(v2, v0)
        l03 = distance(v3, v0) 
        l12 = distance(v2, v1)
        l13 = distance(v3, v1) 
        l23 = distance(v3, v2)
            
        #Volume of cell 
        V = 1/6 * abs(np.linalg.det([e01, e02,e03]))
        D = sqrt((l01**2+ l02**2 +l03**2 + l12**2 +l13**2 + l23**2)**3)
        
        return TetraElementQuality.C_tetra *(V/D)
    
    @staticmethod   
    def skewness(v0: np.ndarray,v1: np.ndarray,v2: np.ndarray,v3: np.ndarray)->float:   
        V  = TetraElementQuality.volume(v0,v1,v2,v3)
        l01 = distance(v1, v0)
        l02 = distance(v2, v0)
        l03 = distance(v3, v0) 
        l12 = distance(v2, v1)
        l13 = distance(v3, v1) 
        l23 = distance(v3, v2)
            
        #Circumradius
        p1 = l01*l23
        p2 = l02*l13
        p3 = l03*l12
        R = sqrt((p1+p2+p3)*(p1-p2+p3)*(p1+p2-p3)*(-p1+p2+p3))/(24*V)
            
        V_ideal = (8*sqrt(3)/27) * (R**3)  
        
        return 1 - V/V_ideal
    
    @staticmethod   
    def edge_ratio(v0: np.ndarray,v1: np.ndarray,v2: np.ndarray,v3: np.ndarray)->float:   
        edge_length = np.array([
            np.linalg.norm(v1-v0),
            np.linalg.norm(v2-v0),
            np.linalg.norm(v3-v0),
            np.linalg.norm(v2-v1),
            np.linalg.norm(v2-v3),
            np.linalg.norm(v3-v1)   
        ])
        
        return max(edge_length)/min(edge_length)
        
    
    #Temp: Needs review and Change!    
    @staticmethod   
    def aspect_ratio_paraview(v0: np.ndarray,v1: np.ndarray,v2: np.ndarray,v3: np.ndarray)->float:   
        edge_length = np.array([
            np.linalg.norm(v1-v0),
            np.linalg.norm(v2-v0),
            np.linalg.norm(v3-v0),
            np.linalg.norm(v2-v1),
            np.linalg.norm(v2-v3),
            np.linalg.norm(v3-v1)   
        ])
        V = TetraElementQuality.volume(v0,v1,v2,v3) 
        
        A = TriElementQuality.area(v0,v1,v2) 
        + TriElementQuality.area(v0,v1,v3) 
        + TriElementQuality.area(v0,v2,v3) 
        + TriElementQuality.area(v1,v2,v3)

        r = 3*V/A
        
        return 4*max(edge_length)/(2*sqrt(6)*r)
      

        
        
        
class TetrahedronMeshQuality():

    C_tetra = 124.70765802
    

    @staticmethod
    def volume(mesh: VolumeMesh)-> np.array:

        volume = np.zeros((mesh.num_cells))
        volume.fill(np.NaN)
        
        volume =np.zeros((mesh.num_cells))
        volume.fill(np.NaN)
    
        for c,cell in enumerate(mesh.cells.astype(int)):
            v0 = mesh.vertices[cell[0]]
            v1 = mesh.vertices[cell[1]]
            v2 = mesh.vertices[cell[2]]
            v3 = mesh.vertices[cell[3]]
            
            volume[c] =  TetraElementQuality.volume(v0,v1,v2,v3)
            
        return volume
    
    @staticmethod
    def element_quality(mesh: VolumeMesh)-> np.array:
        quality = np.zeros((mesh.num_cells))
        quality.fill(np.NaN)
        
        for c,cell in enumerate(mesh.cells.astype(int)):
            v0 = mesh.vertices[cell[0]]
            v1 = mesh.vertices[cell[1]]
            v2 = mesh.vertices[cell[2]]
            v3 = mesh.vertices[cell[3]]
        
            quality[c] = TetraElementQuality.element_quality(v0,v1,v2,v3)
            
        return quality
    

    @staticmethod
    def skewness(mesh: VolumeMesh)-> np.array:
        skewness = np.zeros((mesh.num_cells))
        skewness.fill(np.NaN)
    
    
        for c,cell in enumerate(mesh.cells.astype(int)):
            v0 = mesh.vertices[cell[0]]
            v1 = mesh.vertices[cell[1]]
            v2 = mesh.vertices[cell[2]]
            v3 = mesh.vertices[cell[3]]
        
            skewness[c] = TetraElementQuality.skewness(v0,v1,v2,v3)

        return np.clip(skewness,0,1)
    

    @staticmethod
    def face_skewness(mesh: VolumeMesh)-> np.array:
        face_skewness = []
        for c,cell in enumerate(mesh.cells.astype(int)):
            v0 = mesh.vertices[cell[0]]
            v1 = mesh.vertices[cell[1]]
            v2 = mesh.vertices[cell[2]]
            v3 = mesh.vertices[cell[3]]
        
            face_skewness.append(TriElementQuality.skewness(v0,v1,v2))
            face_skewness.append(TriElementQuality.skewness(v0,v1,v3))
            face_skewness.append(TriElementQuality.skewness(v0,v2,v3))
            face_skewness.append(TriElementQuality.skewness(v1,v2,v3))

        return np.clip(face_skewness,0,1)


    @staticmethod
    def orthogonal_quality(mesh: VolumeMesh)-> np.array:
        neighbors ,shared_faces = get_cell_adj_cells(mesh)
       
        for c,cell in enumerate(mesh.cells.astype(int)):
            v0 = mesh.vertices[cell[0]]
            v1 = mesh.vertices[cell[1]]
            v2 = mesh.vertices[cell[2]]
            v3 = mesh.vertices[cell[3]]

            centroid = (v0+v1+v2+v3)/4
            orthogonality = []
            for v in range(4):
                face = np.delete(cell, v)
                face_centroid = ( mesh.vertices[face[0]]
                                    + mesh.vertices[face[1]]
                                    + mesh.vertices[face[2]]
                                    ) / 3
                A = unit_vector( np.cross(
                    mesh.vertices[face[1]] - mesh.vertices[face[0]],
                    mesh.vertices[face[2]] - mesh.vertices[face[0]]))

                f = unit_vector(face_centroid - centroid)
                
                orthogonality.append(abs(np.dot(A,f)))
                
            for adj_cell,adj_faces in zip(neighbors[c] ,shared_faces[c]):
                
                adj_centroid = (  mesh.vertices[mesh.cells[adj_cell][0]]
                                + mesh.vertices[mesh.cells[adj_cell][1]]
                                + mesh.vertices[mesh.cells[adj_cell][2]]
                                + mesh.vertices[mesh.cells[adj_cell][3]]) / 4

                A = unit_vector( np.cross(
                    mesh.vertices[adj_faces[1]] - mesh.vertices[adj_faces[0]],
                    mesh.vertices[adj_faces[2]] - mesh.vertices[adj_faces[0]]))
            
                c = unit_vector(adj_centroid  - centroid)
                orthogonality.append(abs(np.dot(A,c)))
            
            skewness = TetraElementQuality.skewness(v0,v1,v2,v3)
            
            orthogonality = min(orthogonality)
        
        return min(orthogonality,1-skewness)
              
 
    @staticmethod 
    def face_aspect_ratio(mesh: VolumeMesh)-> np.array:
        face_aspect_ratio = []
        
        for c,cell in enumerate(mesh.cells.astype(int)):
            v0 = mesh.vertices[cell[0]]
            v1 = mesh.vertices[cell[1]]
            v2 = mesh.vertices[cell[2]]
            v3 = mesh.vertices[cell[3]]
        

            face_aspect_ratio.append(TriElementQuality.aspect_ratio(v0,v1,v2))
            face_aspect_ratio.append(TriElementQuality.aspect_ratio(v0,v1,v3))
            face_aspect_ratio.append(TriElementQuality.aspect_ratio(v0,v2,v3))
            face_aspect_ratio.append(TriElementQuality.aspect_ratio(v1,v2,v3))
        return face_aspect_ratio       
    
    @staticmethod     
    def edge_ratio(mesh: VolumeMesh)-> np.array:
        edge_ratio = np.zeros((mesh.num_cells))
        edge_ratio.fill(np.NaN)
        
        for c,cell in enumerate(mesh.cells.astype(int)):
            v0 = mesh.vertices[cell[0]]
            v1 = mesh.vertices[cell[1]]
            v2 = mesh.vertices[cell[2]]
            v3 = mesh.vertices[cell[3]]
        
            edge_ratio[c] = TetraElementQuality.edge_ratio(v0,v1,v2,v3)
            
        return edge_ratio    
    
    
    @staticmethod 
    def face_edge_ratio(mesh: VolumeMesh)-> np.array:
        face_edge_ratio = []
        
        for c,cell in enumerate(mesh.cells.astype(int)):
            v0 = mesh.vertices[cell[0]]
            v1 = mesh.vertices[cell[1]]
            v2 = mesh.vertices[cell[2]]
            v3 = mesh.vertices[cell[3]]
        

            face_edge_ratio.append(TriElementQuality.edge_ratio(v0,v1,v2))
            face_edge_ratio.append(TriElementQuality.edge_ratio(v0,v1,v3))
            face_edge_ratio.append(TriElementQuality.edge_ratio(v0,v2,v3))
            face_edge_ratio.append(TriElementQuality.edge_ratio(v1,v2,v3))
        return face_edge_ratio           


    @staticmethod    
    def size_change(mesh: VolumeMesh)-> np.array:
        neighbors , _ = get_cell_adj_cells(mesh)
        
        size_change = []
        volume = TetrahedronMeshQuality.volume(mesh)
        
        for c in range(mesh.num_cells):
            for adj_cell in neighbors[c]:
                size_change.append(volume[c]/volume[adj_cell])
                
        return size_change
    
    
    @staticmethod   
    def area(mesh: VolumeMesh)-> np.array:
        area = np.zeros((mesh.num_cells,4))
        area.fill(np.NaN)
        
        for c,cell in enumerate(mesh.cells.astype(int)):
            v0 = mesh.vertices[cell[0]]
            v1 = mesh.vertices[cell[1]]
            v2 = mesh.vertices[cell[2]]
            v3 = mesh.vertices[cell[3]]

            area[c] = [TriElementQuality.area(v0,v1,v2),
                       TriElementQuality.area(v0,v1,v3),
                       TriElementQuality.area(v0,v2,v3),
                       TriElementQuality.area(v1,v2,v3)]
            
            
        return area
    
    
def main():
    path_to_mesh = "./sandbox/MeshQualityMetrics/test_meshes/unitCube.inp"
    volume_mesh = dtcc_io.load_volume_mesh(path_to_mesh)
    

    # Single ideal Tetrahedron cell for testing
    # volume_mesh = VolumeMesh()
    # volume_mesh.vertices = np.array([
    # [0.0, 0.0, 0.0],
    # [1.0, 0.0, 0.0],
    # [0.5, np.sqrt(3) / 2, 0.0],
    # [0.5, np.sqrt(3) / 6, np.sqrt(6) / 3]
    # ])
    # volume_mesh.cells = np.array([[0,1,2,3]]) 

    
    a = TetrahedronMeshQuality.area(volume_mesh)
    print("\nFace Area")
    print("MAX\t",       np.nanmax(a))
    print("MIN\t",       np.nanmin(a))
    print("MEAN\t",     np.nanmean(a))
    print("MEDIAN\t", np.nanmedian(a))

    v = TetrahedronMeshQuality.volume(volume_mesh)
    print("\nTetrahedron Element Volume")
    print("MAX\t",       np.nanmax(v))
    print("MIN\t",       np.nanmin(v))
    print("MEAN\t",     np.nanmean(v))
    print("MEDIAN\t", np.nanmedian(v))
     
    eq = TetrahedronMeshQuality.element_quality(volume_mesh)
    print("\nTetrahedron Element Quality")
    print("MAX\t",       np.nanmax(eq))
    print("MIN\t",       np.nanmin(eq))
    print("MEAN\t",     np.nanmean(eq))
    print("MEDIAN\t", np.nanmedian(eq))
    
    sk = TetrahedronMeshQuality.skewness(volume_mesh)
    print("\nOptimal Cell Size Skewness")
    print("MAX\t",       np.nanmax(sk))
    print("MIN\t",       np.nanmin(sk))
    print("MEAN\t",     np.nanmean(sk))
    print("MEDIAN\t", np.nanmedian(sk))
    
    sk = TetrahedronMeshQuality.face_skewness(volume_mesh)
    print("\nOptimal Cell Size Skewness")
    print("MAX\t",       np.nanmax(sk))
    print("MIN\t",       np.nanmin(sk))
    print("MEAN\t",     np.nanmean(sk))
    print("MEDIAN\t", np.nanmedian(sk))
    
    q = TetrahedronMeshQuality.size_change(volume_mesh)
    print("\nCell Volume change")
    print("MAX\t",       np.nanmax(q))
    print("MIN\t",       np.nanmin(q))
    print("MEAN\t",     np.nanmean(q))
    print("MEDIAN\t", np.nanmedian(q))
    
    oq = TetrahedronMeshQuality.orthogonal_quality(volume_mesh)
    print("\nOrthogonality")
    print("MAX\t",       np.nanmax(oq))
    print("MIN\t",       np.nanmin(oq))
    print("MEAN\t",     np.nanmean(oq))
    print("MEDIAN\t", np.nanmedian(oq))

    ar = TetrahedronMeshQuality.face_aspect_ratio(volume_mesh)
    print("\nAspect Ratio (face)")
    print("MAX\t",       np.nanmax(ar))
    print("MIN\t",       np.nanmin(ar))
    print("MEAN\t",     np.nanmean(ar))
    print("MEDIAN\t", np.nanmedian(ar))
    
    
    er = TetrahedronMeshQuality.edge_ratio(volume_mesh)
    print("\nEdge Ratio")
    print("MAX\t",       np.nanmax(er))
    print("MIN\t",       np.nanmin(er))
    print("MEAN\t",     np.nanmean(er))
    print("MEDIAN\t", np.nanmedian(er))
    
    fer = TetrahedronMeshQuality.face_edge_ratio(volume_mesh)
    print("\nFace Edge Ratio")
    print("MAX\t",       np.nanmax(fer))
    print("MIN\t",       np.nanmin(fer))
    print("MEAN\t",     np.nanmean(fer))
    print("MEDIAN\t", np.nanmedian(fer))
    

    

if __name__ == "__main__":
    main()