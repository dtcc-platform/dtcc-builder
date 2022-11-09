from typing import List, Tuple


def bounds_union(
    a: Tuple[float, float, float, float], b: Tuple[float, float, float, float]
):
    px = min(a[0], b[0])
    qx = max(a[2], b[2])
    py = min(a[1], b[1])
    qy = max(a[3], b[3])
    return (px, py, qx, qy)


def bounds_intersect(
    a: Tuple[float, float, float, float], b: Tuple[float, float, float, float]
):
    px = max(a[0], b[0])
    qx = min(a[2], b[2])
    py = max(a[1], b[1])
    qy = min(a[3], b[3])
    return (px, py, qx, qy)


def bounds_area(b: Tuple[float, float, float, float]) -> float:
    return (b[2] - b[0]) * (b[3] - b[1])
