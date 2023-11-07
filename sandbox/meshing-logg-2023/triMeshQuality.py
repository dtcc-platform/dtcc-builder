# Latest Version 21/10/23

import numpy as np
import dtcc_io
import time
from dtcc_builder.logging import *
from dtcc_model import Mesh, VolumeMesh
from qualityMetricsUtils import *

class TriElementQuality():

    C_tri = 6.92820323

    @staticmethod
    def area(v0: np.ndarray,v1: np.ndarray,v2: np.ndarray)->float:
        l0 = np.linalg.norm(v2-v1)
        l1 = np.linalg.norm(v0-v2)
        l2 = np.linalg.norm(v1-v0)

        #Triangle face Semiperimeter
        s = (l0 + l1 + l2)/2
        return sqrt(s*(s-l0)*(s-l1)*(s-l2))

    @staticmethod
    def element_quality(v0: np.ndarray,v1: np.ndarray,v2: np.ndarray)->float:
        l0 = np.linalg.norm(v2-v1)
        l1 = np.linalg.norm(v0-v2)
        l2 = np.linalg.norm(v1-v0)

        area = TriElementQuality.area(v0,v1,v2)
        return TriElementQuality.C_tri * area / (l0**2+l1**2+l2**2)


    @staticmethod
    def edge_ratio(v0: np.ndarray,v1: np.ndarray,v2: np.ndarray)->float:
        l0 = np.linalg.norm(v2-v1)
        l1 = np.linalg.norm(v0-v2)
        l2 = np.linalg.norm(v1-v0)

        l_max = max((l0,l1,l2))
        l_min = min((l0,l1,l2))

        if l_min == 0:
            return np.inf
        return l_max/l_min

    @staticmethod
    def aspect_ratio(v0: np.ndarray,v1: np.ndarray,v2: np.ndarray)->float:
        l0 = np.linalg.norm(v2-v1)
        l1 = np.linalg.norm(v0-v2)
        l2 = np.linalg.norm(v1-v0)

        l_max = max((l0,l1,l2))

        Area = 0.5 * np.linalg.norm(np.cross(v1-v0,v2-v0))

        return l_max*(l0+l1+l2)/(4*sqrt(3)*Area)



    @staticmethod
    def equiangular_skew(v0: np.ndarray,v1: np.ndarray,v2: np.ndarray)->float:

            theta_0 = angle_between(v1-v0,v2-v0)
            theta_1 = angle_between(v2-v1,v0-v1)
            theta_2 = angle_between(v1-v2,v0-v2)

            theta_max = max((theta_0,theta_1,theta_2))
            theta_min = min((theta_0,theta_1,theta_2))

            return max(
            (theta_max - np.pi / 3) / (2 * np.pi / 3),
            (np.pi / 3 - theta_min) / (np.pi / 3)
            )

    @staticmethod
    def skewness(v0: np.ndarray,v1: np.ndarray,v2: np.ndarray)->float:
            l0 = np.linalg.norm(v2-v1)
            l1 = np.linalg.norm(v0-v2)
            l2 = np.linalg.norm(v1-v0)

            #Triangle face Semiperimeter
            s = (l0 + l1 + l2)/2
            area = sqrt(s*(s-l0)*(s-l1)*(s-l2))

            R = (l0*l1*l2) / (4*area)

            area_ideal = (3*sqrt(3)/4)*(R**2)

            skew = 1 - area/area_ideal

            return np.clip(skew,0,1)







# The TriMeshQuality class is used to calculate the quality of a triangular mesh.
class TriMeshQuality():

    C_tri = 6.92820323

    @staticmethod
    def aspect_ratio(mesh: Mesh)->np.ndarray:

        aspect_ratio = np.zeros((mesh.num_faces))
        aspect_ratio.fill(np.NaN)

        for i, face in enumerate(mesh.faces):
            v0 = mesh.vertices[face[0]]
            v1 = mesh.vertices[face[1]]
            v2 = mesh.vertices[face[2]]

            aspect_ratio[i] = TriElementQuality.aspect_ratio(v0,v1,v2)

        return np.around(aspect_ratio, decimals = 8)

    @staticmethod
    def edge_ratio(mesh: Mesh)->np.ndarray:

        edge_ratio = np.zeros((mesh.num_faces))
        edge_ratio.fill(np.NaN)

        for i, face in enumerate(mesh.faces):
            v0 = mesh.vertices[face[0]]
            v1 = mesh.vertices[face[1]]
            v2 = mesh.vertices[face[2]]

            edge_ratio[i] = TriElementQuality.edge_ratio(v0,v1,v2)

        return np.around(edge_ratio, decimals = 8)

    @staticmethod
    def equiangular_skew(mesh: Mesh)->np.ndarray:
        skew = np.zeros((mesh.num_faces))
        skew.fill(np.NaN)

        for i,face in enumerate(mesh.faces):
            v0 = mesh.vertices[face[0]]
            v1 = mesh.vertices[face[1]]
            v2 = mesh.vertices[face[2]]

            skew[i] =  TriElementQuality.equiangular_skew(v0,v1,v2)

        return np.clip(skew,0,1)

    @staticmethod
    def skewness(mesh: Mesh )->np.ndarray:
        skew = np.zeros((mesh.num_faces))
        skew.fill(np.NaN)

        for i,face in enumerate(mesh.faces):
            v0 = mesh.vertices[face[0]]
            v1 = mesh.vertices[face[1]]
            v2 = mesh.vertices[face[2]]

            skew[i] = TriElementQuality.skewness(v0,v1,v2)

        return np.clip(skew,0,1)

    @staticmethod
    def element_quality(mesh: Mesh)->np.ndarray:
        quality = np.zeros((mesh.num_faces), dtype=float)
        quality.fill(np.NaN)

        for i,face in enumerate(mesh.faces):
            v0 = mesh.vertices[face[0]]
            v1 = mesh.vertices[face[1]]
            v2 = mesh.vertices[face[2]]
            quality[i] =  TriElementQuality.element_quality(v0,v1,v2)

        return np.clip(quality,0,1)


    @staticmethod
    def relative_size_squared(mesh: Mesh)-> np.array:
        area = TriMeshQuality.area(mesh)
        mean_area = np.mean(area)
        rss = np.zeros((mesh.num_faces), dtype=float)
        for i,face in enumerate(mesh.faces):
            R = area[i]/mean_area
            rss[i] = min((R,1/R))**2

        return rss

    @staticmethod
    def area(mesh: Mesh)-> np.array:
        area = np.zeros((mesh.num_faces), dtype=float)

        for i,face in enumerate(mesh.faces):
            v0 = mesh.vertices[face[0]]
            v1 = mesh.vertices[face[1]]
            v2 = mesh.vertices[face[2]]

            area[i] = TriElementQuality.area(v0,v1,v2)
        return area

def check_surface_mesh(mesh: Mesh):

    title("Quality metrics for surface mesh")

    a = TriMeshQuality.area(mesh)
    title("Face area (allowed to vary)")
    print("MIN\t",       np.nanmin(a))
    print("MAX\t",       np.nanmax(a))
    print("MEAN\t",     np.nanmean(a))
    print("MEDIAN\t", np.nanmedian(a))

    eq = TriMeshQuality.element_quality(mesh)
    title("Element quality (0-1, 1 is optimal)")
    print("MIN\t",       np.nanmin(eq))
    print("MAX\t",       np.nanmax(eq))
    print("MEAN\t",     np.nanmean(eq))
    print("MEDIAN\t", np.nanmedian(eq))

    er = TriMeshQuality.edge_ratio(mesh)
    title("Edge ratio (1-inf, 1 is optimal)")
    print("MIN\t",       np.nanmin(er))
    print("MAX\t",       np.nanmax(er))
    print("MEAN\t",     np.nanmean(er))
    print("MEDIAN\t", np.nanmedian(er))

    ar = TriMeshQuality.aspect_ratio(mesh)
    title("Aspect ratio (1-inf, 1 is optimal)")
    print("MIN\t",       np.nanmin(ar))
    print("MAX\t",       np.nanmax(ar))
    print("MEAN\t",     np.nanmean(ar))
    print("MEDIAN\t", np.nanmedian(ar))

    sk = TriMeshQuality.equiangular_skew(mesh)
    title("Equiangular skewness (0-1, 0 is optimal)")
    print("MIN\t",       np.nanmin(sk))
    print("MAX\t",       np.nanmax(sk))
    print("MEAN\t",     np.nanmean(sk))
    print("MEDIAN\t", np.nanmedian(sk))

    sk = TriMeshQuality.skewness(mesh)
    title("Skewness (0-1, 0 is optimal)")
    print("MIN\t",       np.nanmin(sk))
    print("MAX\t",       np.nanmax(sk))
    print("MEAN\t",     np.nanmean(sk))
    print("MEDIAN\t", np.nanmedian(sk))

    # This metric is not relevant?
    #rss = TriMeshQuality.relative_size_squared(mesh)
    #title("Relative size squared")
    #print("MIN\t",       np.nanmin(rss))
    #print("MAX\t",       np.nanmax(rss))
    #print("MEAN\t",     np.nanmean(rss))
    #print("MEDIAN\t", np.nanmedian(rss))



def main():

    path_to_mesh = "./surface_mesh.obj"
    mesh = dtcc_io.load_mesh(path_to_mesh)


    # Single Equitriangular Cell for Testing
    # mesh = Mesh()
    # mesh.vertices = np.array([[0, 0, 0], [1, 0, 0], [0.5, np.sqrt(3)/2, 0]])
    # mesh.faces = np.array([[0,1,2]])


    check_surface_mesh(mesh)


if __name__ == "__main__":
    main()