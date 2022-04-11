import numpy as np
import random, math
from scipy.spatial import ConvexHull
import point_cloud_utils as pcu
from sklearn.neighbors import NearestNeighbors
from shapely.geometry import Point, MultiPoint, LineString, Polygon

def get_angle_between_normals(x,y):
    if all([np.float32(xs)==np.float32(ys)  for xs,ys in zip(x,y)]):
        return 0
    else:
        try:
            r = np.dot(x, y) / (np.sqrt(np.sum(np.square(x))) * np.sqrt(np.sum(np.square(y))))
            if r > 1:
                r = 1
            return math.degrees(math.acos(r))
        except BaseException as e:
            print(x,y)
            return 0


def region_growing_segmentation_3d(point_cloud:np.ndarray,max_radius=0.2, normals_angle_threshold = 2):
    point_cloud_normalized = (point_cloud - point_cloud.min(axis=0)) / point_cloud.std(axis=0)

    nbrs = NearestNeighbors(n_neighbors = 20, radius=2)
    nbrs.fit(point_cloud_normalized)
    distances, indexes = nbrs.kneighbors(point_cloud_normalized)

    ## NOTE maybe downsample if needed?
    # pcu.downsample_point_cloud_voxel_grid()

    _,normals = pcu.estimate_point_cloud_normals_knn(np.array(point_cloud,dtype=np.float64),16)

    n_points = len(point_cloud_normalized)
    
    seeds_size = int(n_points*0.1)
    seed_indexes = random.sample(list(np.arange(0,n_points)),seeds_size)
    
    regions_indices = []
    assigned_points = []
    for ls in list(seed_indexes):
        S = set()
        S.add(ls)
        R = []
        while len(S)>0:
            p = S.pop()
            p_normal_vector = normals[p]
            for d,ix in zip(distances[p],indexes[p]):
                if ix not in assigned_points:
                    angle_between_normals = get_angle_between_normals(p_normal_vector,normals[ix])
                    if d<=max_radius and angle_between_normals < normals_angle_threshold:
                        S.add(ix)
                        R.append(int(ix))
                        assigned_points.append(ix)
        if len(R)>0:
            regions_indices.append(R)



    regions = [(point_cloud[region_indices],normals[region_indices],region_indices) for region_indices in regions_indices]
    return regions

def filter_regions(point_cloud, regions, min_points_per_region=10, map_unsegmented =True):
    xy_normal_vector = [0,0,1]
    filtered_regions = [{"size":len(r[0]),"angle":float(get_angle_between_normals(np.mean(r[1],axis=0),xy_normal_vector)), "indices":list(r[2])} for r in regions if len(r[0])>min_points_per_region]

    if map_unsegmented:
        pts = set(range(len(point_cloud)))
        if len(filtered_regions)==0:
            return [{"size":len(pts),"angle":0.0,"indices":sorted(list(pts))}]
        segmented_points = set()
        segments = []
        for region in filtered_regions:
            segmented_points = segmented_points.union(set(region["indices"]))
            segments.append(region["indices"])
        unsegmented_points = pts - segmented_points
        segment_hulls = [MultiPoint(point_cloud[s]).convex_hull for s in segments]
        for pt_idx in unsegmented_points:
            pt = Point(point_cloud[pt_idx])
            closest = np.argmin([pt.distance(s) for s in segment_hulls])
            segments[closest].append(pt_idx)
        for idx,region in enumerate(filtered_regions):
            region["indices"] = sorted(segments[idx])
            region["size"] = len(region["indices"])

    return filtered_regions

def get_roof_segments(point_cloud, min_points_per_region=10, map_unsegmented =True):
    regions = region_growing_segmentation_3d(point_cloud)
    regions = filter_regions(point_cloud,regions,min_points_per_region, map_unsegmented)
    return regions




