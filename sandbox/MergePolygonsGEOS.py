# Prototyping algorithm for merging polygons (union) using GEOS
#
# Note: GEOS = Geometry Engine, Open Source (libgeos), not to be
# confused with GEOS = Google Earth Overlay Server... Name of the
# Python package is pygeos (not geos) so do pip install pygeos.
#
# Anders Logg 2022-04-07

from pylab import *
from pygeos import *

ORANGE = '#F26C24'
BLUE = '#1F95BA'
GREY = '#3a3a3a'

TOL = 0.1
EPS = 1e-06

def Polygon(coordinates):
    return polygons(coordinates)

def TestCase0():
    p0 = array([(0, 0), (1, 0), (1, 1), (0, 1)])
    p1 = p0 + array((1.1, 0.5))
    return polygons(p0), polygons(p1), 'Test case 0'

def TestCase1():
    p0 = array([(0, 0), (1, 0), (1, 1), (0, 1)])
    p1 = p0 + array((0.6, 0.5))
    return polygons(p0), polygons(p1), 'Test case 1'

def TestCase2():
    p0 = array([(0, 0), (1, 0), (1, 1), (0, 1)])
    p1 = array([(1.1, 0.5), (1.5, 0), (2, 0.5), (1.5, 1)])
    return polygons(p0), polygons(p1), 'Test case 2'

def TestCase3():
    p0 = array([(0, 0), (1, 0), (1, 1), (0, 1)])
    p1 = p0 + array((1.0, 0.0))
    return polygons(p0), polygons(p1), 'Test case 3'

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
    return polygons(p0), polygons(p1), 'Test case 4'

def TestCase5():
    p0 = array([(0, 0), (1, 0), (1, 1), (0, 1)])
    p1 = p0 + array((1.1, 0.5))
    p1[0] -= array((0.1, 0))
    return polygons(p0), polygons(p1), 'Test case 5'

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
    return polygons(p0), polygons(p1), 'Test case 6'

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
    return polygons(p0), polygons(p1), 'Test case 7'

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
    return polygons(p0), polygons(p1), 'Test case 8'

def TestCase9():
    p0 = [[539.39399997, 411.2619951 ],
          [534.37199997, 411.5369951 ],
          [533.99599997, 404.7629951 ],
          [539.02399997, 404.4879951 ]]
    p1 = [[526.26499997, 393.8649951 ],
          [553.67099997, 392.2179951 ],
          [554.30599997, 403.8439951 ],
          [533.99599997, 404.7629951 ],
          [534.94199997, 421.9599951 ],
          [523.44399997, 422.4769951 ],
          [522.29099997, 404.43099509],
          [526.81199997, 404.2289951 ]]
    return polygons(p0), polygons(p1), 'Test case 9'

def TestCase10():
    p0 = [[526.26499997, 393.8649951 ],
          [553.67099997, 392.2179951 ],
          [554.30599997, 403.8439951 ],
          [539.02658749, 404.53536783],
          [539.39399997, 411.2619951 ],
          [534.37199997, 411.5369951 ],
          [534.36864511, 411.53717965],
          [534.94199997, 421.9599951 ],
          [523.44399997, 422.4769951 ],
          [522.29099997, 404.43099509],
          [526.81199997, 404.2289951 ]]
    p1 = [[523.16699997, 418.0989951 ],
          [521.35499997, 418.21899509],
          [520.64999997, 407.1569951 ],
          [522.46099997, 407.0369951 ]]
    return polygons(p0), polygons(p1), 'Test case 10'

def TestCase11():
    p0 = [[568.54499997, 407.83799511],
          [554.56599997, 408.6009951 ],
          [552.50099997, 370.9129951 ],
          [566.47899997, 370.14999511]]
    p1 = [[526.26499997, 393.8649951 ],
          [553.67099997, 392.2179951 ],
          [554.30599997, 403.8439951 ],
          [539.02658749, 404.53536783],
          [539.39399997, 411.2619951 ],
          [534.37199997, 411.5369951 ],
          [534.36864511, 411.53717965],
          [534.94199997, 421.9599951 ],
          [523.44399997, 422.4769951 ],
          [523.16429105, 418.09917449],
          [521.35499997, 418.21899509],
          [520.64999997, 407.1569951 ],
          [522.45751798, 407.03722582],
          [522.29099997, 404.43099509],
          [526.81199997, 404.2289951 ]]
    return polygons(p0), polygons(p1), 'Test case 11'

def TestCase12():
    p0 = [[552.50099997, 370.9129951 ],
          [566.47899997, 370.14999511],
          [568.54499997, 407.83799511],
          [554.56599997, 408.6009951 ],
          [554.30535647, 403.84403036]]
    p1 = [[554.42399997, 370.7089951 ],
          [554.40499997, 368.7879951 ],
          [567.26199997, 368.04799511],
          [567.85999997, 380.42599511],
          [567.16399997, 380.36599511],
          [567.03499997, 380.26199511],
          [566.47899997, 370.14499511],
          [554.54399997, 370.8039951 ]]
    return polygons(p0), polygons(p1), 'Test case 12'

def TestCase13():
    p0 = [[292.37499997, 358.47599504],
          [291.75699997, 346.98899504],
          [295.37799997, 346.74599504],
          [295.12199997, 342.06399505],
          [346.97999997, 339.26399506],
          [347.79299997, 354.31499506],
          [330.95999997, 355.22299505],
          [331.11599997, 358.10699505],
          [325.27699997, 358.42299505],
          [325.15126955, 356.09859427],
          [325.12099997, 355.53899505],
          [318.79399997, 355.88299505],
          [308.35499997, 356.44199505],
          [308.42499997, 357.61999505]]
    p1 = [[338.79999997, 325.71299506],
          [338.81999997, 326.14699506],
          [349.68999997, 325.61399506],
          [350.08399997, 329.20399506],
          [350.85899997, 329.17599506],
          [350.91399997, 330.10899506],
          [347.11499997, 334.38699506],
          [347.46599997, 339.26299506],
          [331.63599997, 340.09299505],
          [331.30799997, 333.47099505],
          [334.88199997, 329.60499505],
          [334.81099997, 327.91199505],
          [336.40799997, 325.82299505]]
    return polygons(p0), polygons(p1), 'Test case 13'

def TestCase14():
    p0 = [[231.14399994979613, 150.44799979962409],
          [230.8629999498953, 156.68399979919195],
          [227.78699994989438, 156.25299980025738],
          [228.67399994982406, 150.07099980022758],
          [231.14399994979613, 150.44799979962409]]
    p1 = [[207.7159999498399, 139.2019997993484],
          [217.76999994984362, 140.09099979978055],
          [216.67799994983943, 151.3609997993335],
          [214.94699994981056, 154.25899980030954],
          [215.61399995000102, 155.51899980101734],
          [227.78699994989438, 156.25299980025738],
          [226.24299994989997, 166.89199980068952],
          [207.03799994988367, 164.323999799788],
          [205.54999994981335, 162.44599980022758],
          [207.7159999498399, 139.2019997993484]]
    return polygons(p0), polygons(p1), 'Test case 14'

def GetVertices(geometry):
    vertices = []
    for part in get_parts(geometry):
        points = get_coordinates(part)
        x = [p[0] for p in points]
        y = [p[1] for p in points]
        vertices.append((x, y))
    return vertices

# Not used
def ExplodePolygon(polygon, a):
    c = centroid(polygon)
    r = minimum_bounding_radius(polygon)
    s = 1.0 + a / r
    cx = get_x(c)
    cy = get_y(c)
    points = get_coordinates(polygon)
    x = [cx + s*(p[0] - cx) for p in points]
    y = [cy + s*(p[1] - cy) for p in points]
    return polygons([list(zip(x, y))])[0]

def ComputeVertexProjections(A, B, tol):
    'Project vertices of polygon A on edges of polygon B if close'

    # Get vertices
    vA = get_coordinates(A)
    vB = get_coordinates(B)


    # Create empty list of projections
    projections = []

    # Iterate over vertices in A
    for i in range(len(vA) - 1):

        # Get vertex
        q = vA[i]

        # Iterate over edges in B
        for j in range(len(vB) - 1):

            # Get edge
            p0 = vB[j]
            p1 = vB[j + 1]

            # Project point to edge
            u = q - p0
            v = p1 - p0
            p = p0 + v * dot(u, v) / dot(v, v)

            # Check whether projected point is inside segment. Check either
            # x or y coordinates depending on which is largest (most stable)
            if abs(v[0]) > abs(v[1]):
                inside = min(p0[0], p1[0]) <= p[0] and p[0] <= max(p0[0], p1[0])
            else:
                inside = min(p0[1], p1[1]) <= p[1] and p[1] <= max(p0[1], p1[1])

            # Use either projection or closest vertex
            if inside:
                d = sqrt(dot(q - p, q - p))
                if d < tol:
                    projections.append(p)
                    projections.append(q)
            else:
                d0 = sqrt(dot(q - p0, q - p0))
                d1 = sqrt(dot(q - p1, q - p1))
                if d0 < d1 and d0 < tol:
                    projections.append(p0)
                    projections.append(q)
                elif d1 < tol:
                    projections.append(p1)
                    projections.append(q)

    return projections

def IsValid(C):
    return get_type_id(C) == 3 and minimum_clearance(C) > TOL

def MergePolygons(A, B):

    # Simplify geometries
    A = simplify(A, TOL)
    B = simplify(B, TOL)

    # Compute union
    C = union(A, B, grid_size=EPS)
    C = simplify(C, TOL)

    # Accept if single geometry and good quality
    if IsValid(C):
        print('Using union')
        return C

    # Gradually increase tolerance for merging
    tol = TOL
    maxiter = 3
    for k in range(maxiter):

        # Compute vertex projections
        pAB = ComputeVertexProjections(A, B, tol)
        pBA = ComputeVertexProjections(B, A, tol)
        projections = pAB + pBA
        print(len(projections))

        # Check that we have at least three vertices
        if len(projections) >= 3:

            # Compute convex hull of projections
            P = convex_hull(polygons(projections))

            # Compute union
            C = union(A, B, grid_size=EPS)
            C = union(C, P, grid_size=EPS)
            C = simplify(C, TOL)

            # Accept if single geometry and good quality
            if IsValid(C):
                print('Using extended union')
                return C

        # Increase tolerance
        tol *= 2

    # Try merging convex hulls
    A = convex_hull(A)
    B = convex_hull(B)
    C = union(A, B, grid_size=EPS)
    C = simplify(C, TOL)
    if IsValid(C):
        print('Using union of convex hulls')
        return C

    # Return convex hull
    C = union(A, B)
    C = convex_hull(C)
    C = simplify(C, TOL)

    print('Using convex hull')
    return C

def CheckSum(vm):
    s = 0.0
    for x, y in vm:
        s += sum(x[:-1]) + sum(y[:-1])
    return s

def RunTestCase(testCase):

    # Create polygons
    p0, p1, title = testCase()

    print(title)
    print('-' * len(title))

    # Merge polygons
    pm = MergePolygons(p0, p1)

    # Compute minimum clearance (measure of quality)
    q0 = minimum_clearance(p0)
    q1 = minimum_clearance(p1)
    qm = minimum_clearance(pm)

    print('Minimum clearance p0:', q0)
    print('Minimum clearance p1:', q1)
    print('Minimum clearance pm:', qm)

    # Get vertices (note: multipolygons)
    v0 = GetVertices(p0)
    v1 = GetVertices(p1)
    vm = GetVertices(pm)

    print('Checksum:', CheckSum(vm))
    print(pm)
    print('')

    # Set axis
    xs = v0[0][0] + v1[0][0]
    ys = v0[0][1] + v1[0][1]
    xMin = min(xs)
    xMax = max(xs)
    yMin = min(ys)
    yMax = max(ys)
    dx = 0.1*(xMax - xMin)
    dy = 0.1*(yMax - yMin)
    _axis = [xMin - dx, xMax + dx, yMin - dy, yMax + dy]

    figure()
    for x, y in v0:
        fill(x, y, '-o', alpha=0.75, color=ORANGE)
    for x, y in v1:
        fill(x, y, '-o', alpha=0.75, color=BLUE)
    for x, y in vm:
        plot(x, y, '--o', color=GREY)
    axis(_axis)
    suptitle(title)

def PlotWKT(wkt):
    geometry = Geometry(wkt)
    for x, y in GetVertices(geometry):
        plot(x, y, '-o')

if __name__ == '__main__':

    #RunTestCase(TestCase0)
    #RunTestCase(TestCase1)
    #RunTestCase(TestCase2)
    #RunTestCase(TestCase3)
    #RunTestCase(TestCase4)
    #RunTestCase(TestCase5)
    #RunTestCase(TestCase6)
    #RunTestCase(TestCase7)
    #RunTestCase(TestCase8)
    #RunTestCase(TestCase9)
    #RunTestCase(TestCase10)
    #RunTestCase(TestCase11)
    #RunTestCase(TestCase12)
    #RunTestCase(TestCase13)
    RunTestCase(TestCase14)

    # Debugging of Test case 5
    #PlotWKT('POLYGON ((0 0, 0 1, 1.05 0.995, 1.1 1.5, 2.1 1.5, 2.1 0.5, 1 0.5, 1 0, 0 0))') # Python
    #PlotWKT('POLYGON ((1 1, 1.1 1.5, 2.1 1.5, 2.1 0.5, 1 0.5, 1 0, 0 0, 0 1, 1 1))') # C++

    #PlotWKT('POLYGON ((231.14399994979613 150.44799979962409, 230.8629999498953 156.68399979919195, 227.78699994989438 156.25299980025738, 228.67399994982406 150.07099980022758, 231.14399994979613 150.44799979962409))')
    #PlotWKT('POLYGON ((207.7159999498399 139.2019997993484, 217.76999994984362 140.09099979978055, 216.67799994983943 151.3609997993335, 214.94699994981056 154.25899980030954, 215.61399995000102 155.51899980101734, 227.78699994989438 156.25299980025738, 226.24299994989997 166.89199980068952, 207.03799994988367 164.323999799788, 205.54999994981335 162.44599980022758, 207.7159999498399 139.2019997993484))')

    show()
