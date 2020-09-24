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

def DistancePointPoint(p, q):
    return Norm(p - q)

def NormDot(v, w):
    return Dot(v, w) / (Norm(v)*Norm(w))

def Contains(edge, point, tol):
    p, q = edge
    v = q - p
    if abs(v[0]) > abs(v[1]):
        return min(p[0], q[0]) - tol < point[0] and max(p[0], q[0]) + tol > point[0]
    else:
        return min(p[1], q[1]) - tol < point[1] and max(p[1], q[1]) + tol > point[1]

# FIXME: Can be done much faster using Orient2D
def Intersects(edge0, edge1):
    p = EdgeIntersection(edge0, edge1)
    return Contains(edge0, p, eps) and Contains(edge1, p, eps)

def CheckIntersections(polygon):
    n = len(polygon)
    for i0 in range(n):
        i1 = (i0 + 1) % n
        p0 = polygon[i0]
        p1 = polygon[i1]
        for j0 in range(i0+1, n):
            if min((i0 - j0) % n, (j0 - i0) % n) < 2: continue
            j1 = (j0 + 1) % n
            q0 = polygon[j0]
            q1 = polygon[j1]
            if Intersects((p0, p1), (q0, q1)):
                #print('Bad intersection')
                #figure()
                #plot([p0[0], p1[0]], [p0[1], p1[1]], '-o')
                #plot([q0[0], q1[0]], [q0[1], q1[1]], '-o')
                #show()
                return False
    return True

def CheckDistances(polygon, minimalDistance):
    tol = minimalDistance*minimalDistance
    n = len(polygon)
    for i0 in range(n):
        i1 = (i0 + 1) % n
        p0 = polygon[i0]
        p1 = polygon[i1]
        for j in range(n):
            if min((i0 - j) % n, (j - i0) % n) < 1: continue
            if min((i1 - j) % n, (j - i1) % n) < 1: continue
            q = polygon[j]
            #print(DistanceSegmentPoint(p0, p1, q))
            #figure()
            #PlotPolygons([polygon], style='-')
            #plot([p0[0], p1[0]], [p0[1], p1[1]], '-')
            #plot(q[0], q[1], 'x')
            #show()
            if SquaredDistanceSegmentPoint(p0, p1, q) < tol:
                #print('Bad distance')
                #figure()
                #plot([p0[0], p1[0]], [p0[1], p1[1]], '-o')
                #plot([q0[0], q1[0]], [q0[1], q1[1]], '-o')
                #show()
                return False

    return True

def CheckAngles(polygon, minimalAngle):
    tol = 1.0 - minimalAngle*minimalAngle
    n = len(polygon)
    for i0 in range(n):
        i1 = (i0 + 1) % n
        i2 = (i0 + 2) % n
        p0 = polygon[i0]
        p1 = polygon[i1]
        p2 = polygon[i2]
        u = p1 - p0
        v = p2 - p1
        u2 = Dot(u, u)
        v2 = Dot(v, v)
        dot = Dot(u, v)
        if dot < 0.0 and dot*dot > tol*u2*v2:
            #print('Bad angle')
            return False

    return True

def CheckPolygon(polygon, minimalDistance, minimalAngle):
    if not CheckIntersections(polygon):
        return False
    if not CheckDistances(polygon, minimalDistance):
        return False
    if not CheckAngles(polygon, minimalAngle):
        return False
    return True

def GetPoint(polygon, i):
    return polygon[i % len(polygon)]

def EdgeSign(p0, p1, q):
    v = p1 - p0
    if abs(v[0] > v[1]):
        l = abs(v[0])
        d0 = abs(p0[0] - q[0])
        d1 = abs(p1[0] - q[0])
    else:
        l = abs(v[1])
        d0 = abs(p0[1] - q[1])
        d1 = abs(p1[1] - q[1])
    if d0 > l - eps and d0 > d1:
        return 1
    elif d1 > l - eps and d1 > d0:
        return -1
    else:
        return 0

def Orient2D(p0, p1, q):
    u = p1 - p0
    v = q - p0
    return u[0]*v[1] - u[1]*v[0]

def ConvexHull(points):

    # The convex hull is computed by doing a Graham scan: select an
    # extreme base point, sort remaining points by angle and then
    # add points that create a left turn around the perimeter.

    # Find point with smallest y-coordinate. If y-coordinate is
    # the same, sort by smallest x-coordinate.
    xMin = points[0][0];
    yMin = points[0][1];
    iMin = 0;
    numPoints = len(points)
    for i in range(1, numPoints):
        x = points[i][0];
        y = points[i][1]
        if (y < yMin or (y == yMin and x < xMin)):
            xMin = x;
            yMin = y;
            iMin = i;

    # Set base point
    baseIndex = iMin
    basePoint = points[baseIndex]

    # Compute angles and distances relative to base point
    angles = []
    for i in range(numPoints):

        # Skip base point
        if i == baseIndex:
            continue;

        # Compute angle (negative cosine) and distance
        p = points[i]
        v = p - basePoint
        distance = Norm(v)
        angle = -v[0] / distance if distance > eps else 0.0

        # Store angle and distance along with index (for sorting)
        angles.append((angle, distance, i))

    # Sort by angles (primary) and distance (secondary) to base point
    angles = sorted(angles)

    # Filter out points with unique angles, keeping only furthest point
    filteredIndices = []
    lastAngle = 2.0 # no angle has this value
    for i in range(numPoints - 1):

        # Get data for current point
        currentAngle = angles[i][0]
        currentIndex = angles[i][2]

        # Add point or replace last point
        if abs(currentAngle - lastAngle) > eps:
            filteredIndices.append(currentIndex)
        else:
            filteredIndices[len(filteredIndices) - 1] = currentIndex

        # Update last index
        lastAngle = currentAngle

    # Create stack of points and add first three candidates
    convexHull = []
    convexHull.append(baseIndex)
    convexHull.append(filteredIndices[0])
    convexHull.append(filteredIndices[1])

    # Graham-Scan: Push candidates to stack and pop until
    # we have a left turn
    for i in range(2, len(filteredIndices)):

        # Get next point
        i2 = filteredIndices[i]
        p2 = points[i2]

        # Keep popping from stack until we see a left turn
        while True:

            # Get last two points from stack
            i1 = convexHull.pop()
            i0 = convexHull[-1]
            p0 = points[i0]
            p1 = points[i1]

            # Check orientation, keep p1 if orientation is positive
            if Orient2D(p0, p1, p2) > eps:
                convexHull.append(i1)
                break

        # Push next candidate to stack
        convexHull.append(i2)

    # Extract polygon points from stack
    polygon = []
    while len(convexHull) > 0:
        polygon.append(points[convexHull.pop()])

    # Reverse polygon to make it counter-clockwise
    polygon.reverse()

    return polygon;

def ConnectVertexEdge(i, j0, j1, points, edges, tol, plotting):

    # Get points
    p = points[i]
    q0 = points[j0]
    q1 = points[j1]

    # Connect vertices if close (create new edge)
    connected = False
    if DistancePointPoint(p, q0) < tol:
        edges[i].append(j0)
        edges[j0].append(i)
        connected = True
    if DistancePointPoint(p, q1) < tol:
        edges[i].append(j1)
        edges[j1].append(i)
        connected = True

    # Don't connect vertex to edge if already connected
    if connected:
        return

    # Don't connect vertex to edge if zero length
    v = q1 - q0
    vNorm = Norm(v)
    if vNorm < eps:
        return

    # Connect vertex to edge if close (project)
    if DistanceSegmentPoint(q0, q1, p) < tol:
        v /= vNorm
        r = q0 + Dot(p - q0, v)*v
        k = len(points)
        points.append(r)
        edges.append([i, j0, j1])
        edges[i].append(k)
        edges[j0].append(k)
        edges[j1].append(k)
        if plotting: plot(r[0], r[1], 'x')

def ConnectEdgeEdge(i0, i1, j0, j1, points, edges, tol, plotting):

    # Get points
    p0 = points[i0]
    p1 = points[i1]
    q0 = points[j0]
    q1 = points[j1]

    # Don't look for intersection if almost parallel
    u = p1 - p0
    v = q1 - q0
    if abs(NormDot(u, v)) > 1.0 - eps:
        return

    # Compute edge-edge intersection
    e = (p0, p1)
    f = (q0, q1)
    r = EdgeIntersection(e, f)

    # Connect edges to intersection if close
    if Contains(e, r, tol) and Contains(f, r, tol):
        k = len(points)
        points.append(r)
        sp = EdgeSign(p0, p1, r)
        sq = EdgeSign(q0, q1, r)
        kEdges = []
        if sp == -1 or sp == 0:
            edges[i0].append(k)
            kEdges.append(i0)
        if sp == 0 or sp == 1:
            edges[i1].append(k)
            kEdges.append(i1)
        if sq == -1 or sq == 0:
            edges[j0].append(k)
            kEdges.append(j0)
        if sq == 0 or sq == 1:
            edges[j1].append(k)
            kEdges.append(j1)
        edges.append(kEdges)
        if plotting: plot(r[0], r[1], 'x')

def MergePolygons(polygons, tol=0.5, plotting=False):
    firstPolygon = polygons[0]
    secondPolygon = polygons[1]

    # Get number of points
    m = len(firstPolygon)
    n = len(secondPolygon)
    #print(m, n)

    # Create list of points
    points = [p for p in firstPolygon] + [p for p in secondPolygon]

    # Create directed graph of edges
    edges = []
    for i in range(m):
        edges.append([(i+1) % m])
    for i in range(n):
        edges.append([(i+1) % n + m])

    # Find all pairwise connections between
    # edge i = (i0, i1) and edge j = (j0, j1)
    for i0 in range(m):
        i1 = edges[i0][0]
        for j0 in range(m, m + n):
            j1 = edges[j0][0]

            # Find vertex-edge connections
            ConnectVertexEdge(i0, j0, j1, points, edges, tol, plotting)
            ConnectVertexEdge(i1, j0, j1, points, edges, tol, plotting)
            ConnectVertexEdge(j0, i0, i1, points, edges, tol, plotting)
            ConnectVertexEdge(j1, i0, i1, points, edges, tol, plotting)

            # Find edge-edge connections
            ConnectEdgeEdge(i0, i1, j0, j1, points, edges, tol, plotting)

    # Remove duplicate vertices
    numPoints = len(points)
    vertexMap = [i for i in range(numPoints)]
    removed = [False for i in range(numPoints)]
    for i in range(numPoints):
        for j in range(i+1, numPoints):
            if removed[i]: continue
            if DistancePointPoint(points[i], points[j]) < eps:
                #print('Merging:', i, j)
                edges[i] = edges[i] + edges[j]
                edges[j] = []
                vertexMap[j] = i
                removed[j] = True

    # Replace removed vertices in graph
    for i in range(numPoints):
        for j in range(len(edges[i])):
            edges[i][j] = vertexMap[edges[i][j]]

    # Remove duplicate edges in graph
    for i in range(len(edges)):
        newEdge = []
        for j in edges[i]:
            if j not in newEdge and i != j: newEdge.append(j)
        edges[i] = newEdge

    #for i, e in enumerate(edges):
    #    print(i, e)

    # Write point labels (and make sure they don't overlap)
    if plotting:
        H = 0.0075*max([Norm(p-q) for p in points for q in points])
        for i, p in enumerate(points):
            h = H
            for j, q in enumerate(points[:i]):
                if Norm(p - q) < eps:
                    h += 5*H
            text(p[0] + h, p[1] + H, str(i), va='bottom', ha='left')

    # Find first vertex by looking for an original edge that is to the
    # "right" of all points
    for i in range(m + n):

        # Skip if no outgoing edges
        if len(edges[i]) == 0: continue

        # Get the edge
        j = edges[i][0]
        u = points[j] - points[i]
        u /= Norm(u)

        # Check all points
        ok = True
        for k in range(numPoints):

            # Skip if removed
            if removed[k]: continue

            # Skip if on edge
            if k == i or k == j: continue

            # Check sin of angle (cross product)
            v = points[k] - points[i]
            v /= Norm(v)
            sin = u[0]*v[1] - u[1]*v[0]
            if sin < -eps:
                ok = False
                break

        # Found first edge
        if ok:
            firstVertex = i
            nextVertex = j
            break

    # Keep track of visited vertices
    visited = [False for i in range(len(points))]
    visited[firstVertex] = True
    visited[nextVertex] = True

    # Initialize polygon
    vertices = [firstVertex, nextVertex]

    # Maximum number of step before failure
    maxNumSteps = 2*numPoints

    # Walk graph to build polygon counter-clockwise by picking
    # the right-most turn at each intersection
    for step in range(maxNumSteps):

        #print('Vertices:', vertices)

        # Get previous and current vertex
        i = len(vertices) - 1
        previousVertex = vertices[i - 1]
        currentVertex = vertices[i]

        # Get current edge(s)
        edge = edges[currentVertex]

        #print('Vertex %d: %s' % (currentVertex, str(edges[currentVertex])))

        # Find next vertex
        if len(edge) == 1:

            # If we only have one edge then follow it
            nextVertex = edge[0]

        else:

            # Get previous edge
            u = points[currentVertex] - points[previousVertex]
            d = Norm(u)
            u = u / d

            # Compute angles and distances for outgoing edges
            angles = []
            for k in edge:

                # Skip if already visited (if not first vertex)
                if k != firstVertex and visited[k]: continue

                # Skip if candidate edge intersects previous edges
                intersects = False
                e = (points[currentVertex], points[k])
                for l in range(1, i - 1):
                    f = (points[vertices[l]], points[vertices[l+1]])
                    if Intersects(e, f):
                        intersects = True
                        break
                if intersects: continue

                # Get new edge
                v = points[k] - points[currentVertex]
                d = Norm(v)
                v = v / d

                # Replace actual angle by cheaper but strictly increasing
                # function to avoid needing to call acos() or asin().
                sin = u[0]*v[1] - u[1]*v[0]
                cos = u[0]*v[0] + u[1]*v[1]
                if cos < -1.0 + eps: continue
                a = sin if cos >= 0.0 else (2.0-sin if sin > 0.0 else sin-2.0)
                angles.append((k, a, d))

            # If we have no more vertices to visit, take a step back
            if len(angles) == 0:
                print('No more vertices to visit, stepping back')
                del vertices[i]
                continue

            # Print angles
            #for angle in angles:
            #    print(angle)

            # Find smallest (right-most) angle. First priority is the angle
            # and second priority is the distance (pick closest). Note that
            # we add a small tolerance to ensure we get the closest vertex
            # if the vertices are on the same line.
            minAngle = angles[0]
            for angle in angles[1:]:
                if (angle[1] < minAngle[1] - eps) or \
                   (angle[1] < minAngle[1] + eps and angle[2] < minAngle[2]):
                    minAngle = angle

            # Pick next vertex
            nextVertex = minAngle[0]

        #print(currentVertex, "-->", nextVertex)
        #print('')

        # We are done if we return to the first vertex
        if nextVertex == firstVertex:
            #print('Back to first vertex')
            break

        # Add next vertex
        vertices.append(nextVertex)
        visited[nextVertex] = True

    # If merge failed, return convex hull
    if nextVertex != firstVertex:
        print('Merge failed, falling back to convex hull')
        points = [p for p in firstPolygon] + [p for p in secondPolygon]
        return ConvexHull(points)

    # Extract polygon points
    polygon = [points[i] for i in vertices]

    return array(polygon)

def PlotLabel(polygon, labelText):
    x = [x[0] for x in polygon]
    y = [x[1] for x in polygon]
    text(mean(x), mean(y), labelText, va='center', ha='center')

def PlotPolygons(polygons, style='-o', arrows=False, labels=False):
    for i, polygon in enumerate(polygons):
        x = [x[0] for x in polygon]
        y = [x[1] for x in polygon]
        if labels: PlotLabel(polygon, str(i))
        x = x + [x[0]]
        y = y + [y[0]]
        if arrows:
            for j in range(len(x) - 1):
                Arrow(x[j], y[j], x[j+1], y[j+1])
        plot(x, y, style)
        axis('off')

def RunTestCase(polygons):
    polygons = [array(polygons[0]).astype(float),
                array(polygons[1]).astype(float)]
    figure()
    subplot(2, 1, 1)
    PlotPolygons(polygons)
    polygon = MergePolygons(polygons, plotting=True)
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
    p0 = [[551.02099997, 57.5619951],
          [557.94999997, 184.41399511],
          [545.64399997, 185.08599511],
          [539.38399997, 70.5609951],
          [530.07099997, 71.0789951],
          [529.47099997, 58.7409951]]
    p1 = [[529.47099997, 58.7409951],
          [530.07099997, 71.0789951],
          [460.38899997, 74.87799509],
          [462.36899997, 111.58199508],
          [449.93899997, 112.26099508],
          [447.26099997, 63.22899508]]
    return p0, p1

def TestCase5():
    p0 = array([(0, 0), (1, 0), (1, 1), (0, 1)])
    p1 = p0 + array((1.1, 0.5))
    p1[0] -= array((0.1, 0))
    return p0, p1

def TestCase6():
    p0 = [[338.82099997, 326.20099506],
          [335.96199997, 268.23299506],
          [346.62299997, 267.70599506],
          [346.53899997, 263.96099506],
          [354.08399997, 263.68399506],
          [357.65099997, 328.92299506],
          [350.08599997, 329.20599506],
          [349.69199997, 325.66799506]]
    p1 = [[291.75699997, 346.98899504],
          [295.37799997, 346.74599504],
          [295.12199997, 342.06399505],
          [331.63597417, 340.09247421],
          [331.30799997, 333.47099505],
          [334.88199997, 329.60499505],
          [334.81099997, 327.91199505],
          [336.40799997, 325.82299505],
          [338.79999997, 325.71299506],
          [338.81999997, 326.14699506],
          [349.68999997, 325.61399506],
          [350.08399997, 329.20399506],
          [350.85899997, 329.17599506],
          [350.91399997, 330.10899506],
          [347.11499997, 334.38699506],
          [347.46419003, 339.23785189],
          [347.46599997, 339.26299506],
          [346.98131866, 339.28840791],
          [347.79299997, 354.31499506],
          [330.95999188, 355.22284556],
          [330.95999997, 355.22299505],
          [331.11599997, 358.10699505],
          [325.27699997, 358.42299505],
          [325.23293055, 357.60827579],
          [318.91799997, 357.95799505],
          [318.86933569, 357.14365328],
          [318.86299997, 357.14399505],
          [308.43015642, 357.70677084],
          [308.42599997, 357.70699505],
          [308.42132643 ,357.62019097],
          [292.37499997, 358.47599504]]
    return p0, p1

def TestCase7():
    p0 = [[689.40299996, 59.61399514],
          [689.35699996, 64.63399514],
          [676.57299996, 65.25199513],
          [676.74399996, 59.23399514]]
    p1 = [[676.56899996, 65.37299513],
          [676.57299996, 65.25199513],
          [689.35699996, 64.63399514],
          [689.35399996, 64.92899514],
          [689.96199996, 69.64799514],
          [677.26999996, 71.23799513]]
    return p0, p1

def TestCase8():
    p0 = [[83.58099997, 249.880995],
          [88.10899997, 249.639995],
          [88.28699997, 252.991995],
          [83.75999997, 253.231995]]
    p1 = [[71.24799997, 250.721995],
          [71.48399997, 254.38399499],
          [75.36099997, 254.137995],
          [75.12499997, 250.480995],
          [83.93599997, 256.027995],
          [64.79599997, 257.07299499],
          [64.44799997, 250.89299499],
          [83.59299997, 249.847995]]
    return p0, p1

if __name__ == '__main__':
    #RunTestCase(TestCase0())
    #RunTestCase(TestCase1())
    #RunTestCase(TestCase2())
    #RunTestCase(TestCase3())
    #RunTestCase(TestCase4())
    #RunTestCase(TestCase5())
    #RunTestCase(TestCase6())
    #RunTestCase(TestCase7())
    RunTestCase(TestCase8())
    show()
