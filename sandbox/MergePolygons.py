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

def QuadrantAnglePointPoint(p, q):
    if p[0] > q[0]:
        if p[1] > q[1]:
            return 0
        else:
            return 3
    else:
        if p[1] > q[1]:
            return 1
        else:
            return 2

def QuadrantAnglePointPolygon(p, polygon):

    # Compute angle to first vertex
    q0 = polygon[0]
    v0 = QuadrantAnglePointPoint(q0, p)

    # Sum up total angle
    totalAngle = 0
    for i in range(1, len(polygon) + 1):

      # Compute angle increment
      q1 = polygon[i % len(polygon)]
      v1 = QuadrantAnglePointPoint(q1, p)
      dv = v1 - v0

      # Adjust angle increment for wrap-around
      if dv == 3:
        dv = -1
      elif dv == -3:
        dv = 1
      elif dv == 2 or dv == -2:
        xx = q1[0] - ((q1[1] - p[1]) * ((q0[0] - q1[0]) / (q0[1] - q1[1])))
        if xx > p[0]:
          dv = -dv;

      # Add to total angle and update
      totalAngle += dv;
      q0 = q1;
      v0 = v1;

    return totalAngle;

def PolygonContainsPoint(polygon, p):
    return QuadrantAnglePointPolygon(p, polygon) != 0

def SquaredDistanceSegmentPoint(p0, p1, q):

    # Project point to line
    u = q - p0
    v = p1 - p0
    p = p0 + v * Dot(u, v) / Dot(v, v)

    # Check whether projected point is inside segment. Check either
    # x or y coordinates depending on which is largest (most stable)
    if abs(v[0]) > abs(v[1]):
        inside = min(p0[0], p1[0]) <= p[0] and p[0] <= max(p0[0], p1[0])
    else:
        inside = min(p0[1], p1[1]) <= p[1] and p[1] <= max(p0[1], p1[1])

    # Use distance to projection if inside
    if inside:
      return Dot(q - p, q - p)

    # Otherwise use distance to closest end point
    d0 = Dot(q - p0, q - p0)
    d1 = Dot(q - p1, q - p1)
    return min(d0, d1)

def DistanceSegmentPoint(p0, p1, q):
    return sqrt(SquaredDistanceSegmentPoint(p0, p1, q))

def SquaredDistancePolygonPoint(polygon, p):

    # Check if point is contained in polygon
    if PolygonContainsPoint(polygon, p):
      return 0.0

    # If not, compute minimal squared distance to all segments
    d2Min = sys.float_info.max
    for i in range(len(polygon)):
        p0 = polygon[i]
        p1 = polygon[(i + 1) % len(polygon)]
        d2Min = min(d2Min, SquaredDistanceSegmentPoint(p0, p1, p))

    return d2Min

def SquaredDistancePolygonPolygon(polygon0, polygon1):

    d2Min = sys.float_info.max

    # Check all vertices in first polygon
    for p in polygon0:
      d2Min = min(d2Min, SquaredDistancePolygonPoint(polygon1, p))

    # Check all vertices in second polygon
    for p in polygon1:
      d2Min = min(d2Min, SquaredDistancePolygonPoint(polygon0, p))

    return d2Min

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

def MergePolygons(polygons, tol=0.2):
    polygon0 = polygons[0]
    polygon1 = polygons[1]

    # Get number of points
    n0 = len(polygon0)
    n1 = len(polygon1)

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

            # First check if a vertex is incident on an edge
            incident = False
            if DistanceSegmentPoint(p1, q1, p0) < eps:
                edges[i0].append(i1)
                edges[i0].append(j1)
                edges[i1].append(i0)
                edges[j1].append(i0)
                incident = True
            if DistanceSegmentPoint(p1, q1, q0) < eps:
                edges[j0].append(i1)
                edges[j0].append(j1)
                edges[i1].append(j0)
                edges[j1].append(j0)
                incident = True
            if DistanceSegmentPoint(p0, q0, p1) < eps:
                edges[i1].append(i0)
                edges[i1].append(j0)
                edges[i0].append(i1)
                edges[j0].append(i1)
                incident = True
            if DistanceSegmentPoint(p0, q0, q1) < eps:
                edges[j1].append(i0)
                edges[j1].append(j0)
                edges[i0].append(j1)
                edges[j0].append(j1)
                incident = True

            # Don't look for intersection if incident
            if incident: continue

            # Don't look for intersection if almost parallel
            if abs(NormDot(v0, v1)) > 0.9: continue

            # Compute intersection of lines defined by edges
            p = EdgeIntersection(e0, e1)
            k = len(points)
            e = []

            # Skip intersection if not close to both edges
            if not (Contains(e0, p, tol) and Contains(e1, p, tol)):
                continue

            # Check first edge
            if Contains(e0, p, eps):
                edges[i0].append(k)
                edges[j0].append(k)
                e.append(i0)
                e.append(j0)
            else:
                di = Norm(p0 - p)
                dj = Norm(q0 - p)
                if di < dj:
                    edges[i0].append(k)
                    e.append(i0)
                else:
                    edges[j0].append(k)
                    e.append(j0)

            # Check second edge
            if Contains(e1, p, eps):
                edges[i1].append(k)
                edges[j1].append(k)
                e.append(i1)
                e.append(j1)
            else:
                di = Norm(p1 - p)
                dj = Norm(q1 - p)
                if di < dj:
                    edges[i1].append(k)
                    e.append(i1)
                else:
                    edges[j1].append(k)
                    e.append(j1)

            # Add new edge to graph
            edges.append(e)
            points.append(p)

            plot(p[0], p[1], 'x')

    print(edges)

    # Write point labels (and make sure they don't overlap)
    H = 0.02
    for i, p in enumerate(points):
        h = H
        for j, q in enumerate(points[:i]):
            if Norm(p - q) < eps:
                h += 5*H
        text(p[0] + h, p[1] + H, str(i), va='bottom', ha='left')

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
        if p[0] < xMin or (p[0] < xMin + eps and p[1] < yMin):
            xMin = p[0]
            yMin = p[1]
            firstVertex = i
    print('firstVertex = ', firstVertex)

    # FIXME: Testing

    # Check that first vertex has only one outgoing edge
#    if len(edges[firstVertex]) != 1:
#        raise RuntimeError('First vertex should have a single edge')

    # Keep track of visited vertices
    visited = [False for i in range(len(points))]

    # Walk along first edge
    i = firstVertex # current vertex
    j = edges[i][0] # next vertex
    print(i, "-->", j)
    u = points[j] - points[i]
    d = Norm(u)
    u = u / d
    polygon = [points[i]]
    visited[i] = True
    i = j

    # Walk graph to build polygon counter-clockwise by picking
    # the right-most turn at each intersection
    while (True):

        # Add point to polygon if not visited before. Note that this if-case
        # handles the case when a single vertex is close to an edge and we
        # end up visiting the same vertex twice. We then just skip that
        # vertex and move on to the next. Seems to work well.
        if not visited[i]:
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
                if d < eps: continue
                v = v / d
                sin = u[0]*v[1] - u[1]*v[0]
                cos = u[0]*v[0] + u[1]*v[1]
                a = sin if cos >= 0.0 else (2.0-sin if sin > 0.0 else sin-2.0)
                angles.append((vertex, a, d))

            # We are done if we run out of vertices to visit
            if len(angles) == 0:
                break

            print('')
            print('Vertex ', i)
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

        # We are done if we return to the first vertex
        if j == firstVertex:
            break

        # We should not return to the same vertex
        #if visited[i]:
        #    print('Error, already visited vertex')
        #    break

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

def TestCase3():
    p0 = array([(0, 0), (1, 0), (1, 1), (0, 1)])
    p1 = p0 + array((1.0, 0.0))
    return p0, p1

def TestCase4():
    p0 = array([[551.02099997, 57.5619951],
                [557.94999997, 184.41399511],
                [545.64399997, 185.08599511],
                [539.38399997, 70.5609951],
                [530.07099997, 71.0789951],
                [529.47099997, 58.7409951]])
    p1 = array([[529.47099997, 58.7409951],
                [530.07099997, 71.0789951],
                [460.38899997, 74.87799509],
                [462.36899997, 111.58199508],
                [449.93899997, 112.26099508],
                [447.26099997, 63.22899508]])
    return p0, p1

if __name__ == '__main__':
    RunTestCase(TestCase0())
    RunTestCase(TestCase1())
    RunTestCase(TestCase2())
    RunTestCase(TestCase3())
    RunTestCase(TestCase4())
    show()
