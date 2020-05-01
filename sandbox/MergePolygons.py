# Prototyping algorithm for merging polygons (union)
# Anders Logg 2020-04-29

import sys
from pylab import *

eps = 1e-6

def Arrow(x0, y0, x1, y1, color='grey', rad=0.2, arrowstyle='->', size=10):
    annotate('',
             xy=(x1, y1), xycoords='data',
             xytext=(x0, y0), textcoords='data',
             size=size, va="center", ha="center",
             arrowprops=dict(arrowstyle=arrowstyle, color=color,
                             connectionstyle="arc3,rad=%g" % rad))

def EdgeIntersection(edge0, edge1):
    p0, q0 = edge0
    p1, q1 = edge1

    u = p1 - p0
    v = q0 - p0
    w = q1 - p1

    a = v[0]
    b = -w[0]
    c = v[1]
    d = -w[1]
    e = u[0]
    f = u[1]

    det = a * d - b * c
    k = (d * e - b * f) / det
    p = p0 + v * k;

    return p

def Dot(v, w):
    return v[0]*w[0] + v[1]*w[1]

def Norm(v):
    return sqrt(v[0]**2 + v[1]**2)

def NormDot(v, w):
    return Dot(v, w) / (Norm(v)*Norm(w))

def Contains(edge, point, tol):
    p, q = edge
    v = q - p
    if abs(v[0]) > abs(v[1]):
        return min(p[0], q[0]) - tol < point[0] and max(p[0], q[0]) + tol > point[0]
    else:
        return min(p[1], q[1]) - tol < point[1] and max(p[1], q[1]) + tol > point[1]

def GetPoint(polygon, i):
    return polygon[i % len(polygon)]

def MergePolygons(polygons):
    polygon0 = polygons[0]
    polygon1 = polygons[1]

    # Get number of points
    n0 = len(polygon0)
    n1 = len(polygon1)

    tol = 0.25

    # Create list of points
    points = [p for p in polygon0] + [p for p in polygon1]

    # Create directed graph of edges
    edges = []
    for i in range(n0):
        edges.append([(i+1) % n0])
    for i in range(n1):
        edges.append([(i+1) % n1 + n0])

    print(edges)

    # Compute pairwise edge intersections
    intersections = []
    for i0 in range(n0):
        j0 = edges[i0][0]
        p0 = points[i0]
        q0 = points[j0]
        e0 = (p0, q0)
        v0 = q0 - p0
        for i1 in range(n1, n0 + n1):
            j1 = edges[i1][0]
            p1 = points[i1]
            q1 = points[j1]
            e1 = (p1, q1)
            v1 = q1 - p1
            if abs(NormDot(v0, v1)) < 0.5:
                p = EdgeIntersection(e0, e1)
                if Contains(e0, p, tol) and Contains(e1, p, tol):
                    k = len(points)
                    edges[i0].append(k)
                    edges[j0].append(k)
                    edges[i1].append(k)
                    edges[j1].append(k)
                    edges.append([i0, j0, i1, j1])
                    points.append(p)
                    plot(p[0], p[1], 'x')

    print(edges)

    # Plot all points
    h = 0.02
    for i, p in enumerate(points):
        text(p[0] + h, p[1] + h, str(i), va='bottom', ha='left')

    # # Sort edges by distance to first point
    # for i in range(n0+1):
    #     if len(edges[0][i]) > 2:
    #         p0 = edges[0][i][0]
    #         edges[0][i].sort(key=lambda p, p0=p0 : Norm(p - p0))

    # Find starting point = lower left corner
    xMin = sys.float_info.max
    yMin = sys.float_info.max
    firstVertex = 0
    for i, p in enumerate(points):
        if p[0] < xMin and p[1] < yMin:
            xMin = p[0]
            yMin = p[1]
            firstVertex = i

    # Check that first vertex has only one outgoing edge
    if len(edges[firstVertex]) != 1:
        raise RuntimeError('First vertex should have a single edge')

    # Keep track of visited vertices
    visited = [False for i in range(len(points))]

    # Walk along first edge
    i = firstVertex # current vertex
    j = edges[i][0] # next vertex
    u = points[j] - points[i]
    d = Norm(u)
    u = u / d
    polygon = [points[i]]
    visited[i] = True
    i = j

    # Walk graph to build polygon counter-clockwise by picking
    # the right-most turn at each intersection
    while (True):

        # Add point to polygon
        polygon.append(points[i])

        # Get current edge(s)
        edge = edges[i]

        # Decide next vertex
        if len(edge) == 1:

            # If only one edge, just follow
            j = edge[0]

        else:

            # If multiple edges, compute angles and distances.
            # Replace actual angle by cheaper but strictly increasing
            # function to avoid needing to call acos() or asin().
            angles = []
            for vertex in edge:
                if visited[vertex]: continue
                v = points[vertex] - points[i]
                d = Norm(v)
                v = v / d
                sin = u[0]*v[1] - u[1]*v[0]
                cos = u[0]*v[0] + u[1]*v[1]
                a = sin if cos >= 0.0 else (2.0-sin if sin > 0.0 else sin-2.0)
                angles.append((vertex, a, d))

            # Done if we run out of vertices to visit
            if len(angles) == 0:
                break

            print('')
            for angle in angles:
                print(angle)

            # Find smallest (right-most) angle. First priority is the angle
            # and second priority is the distance (pick closest). Note that
            # we add a small tolerance to ensure we get the closest vertex
            # if the vertices are on the same line.
            minAngle = angles[0]
            for angle in angles[1:]:
                if (angle[1] < minAngle[1] - eps) or (angle[1] < minAngle[1] + eps and angle[2] < minAngle[2]):
                    minAngle = angle

            # Pick next vertex
            j = minAngle[0]

        print(i, "-->", j)

        # Com

        # Move to next vertex
        u = points[j] - points[i]
        d = Norm(u)
        u = u / d
        visited[i] = True
        i = j

    return polygon

def PlotPolygons(polygons, style='-o', arrows=False):
    for polygon in polygons:
        x = [x[0] for x in polygon]
        y = [x[1] for x in polygon]
        x = x + [x[0]]
        y = y + [y[0]]
        if arrows:
            for i in range(len(x) - 1):
                Arrow(x[i], y[i], x[i+1], y[i+1])
        plot(x, y, style)
        axis('off')

def RunTestCase(polygons):
    figure()
    subplot(2, 1, 1)
    PlotPolygons(polygons)
    polygon = MergePolygons(polygons)
    subplot(2, 1, 2)
    PlotPolygons([polygon], '--o', arrows=True)

def TestCase0():
    p0 = array([(0, 0), (1, 0), (1, 1), (0, 1)])
    p1 = p0 + array((1.1, 0.5))
    return p0, p1

def TestCase1():
    p0 = array([(0, 0), (1, 0), (1, 1), (0, 1)])
    p1 = p0 + array((0.6, 0.5))
    return p0, p1

def TestCase2():
    p0 = array([(0, 0), (1, 0), (1, 1), (0, 1)])
    p1 = array([(1.1, 0.5), (1.5, 0), (2, 0.5), (1.5, 1)])
    return p0, p1

RunTestCase(TestCase0())
RunTestCase(TestCase1())
#RunTestCase(TestCase2())

show()
