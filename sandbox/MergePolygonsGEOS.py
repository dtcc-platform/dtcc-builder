# Prototyping algorithm for merging polygons (union) using GEOS
#
# Note: GEOS = Geometry Engine, Open Source (libgeos), not to be
# confused with GEOS = Google Earth Overlay Server... Name of the
# Python package is pygeos (not geos) so do pip install pygeos.
#
# Anders Logg 2022-04-07

from pylab import *
from pygeos import *

def TestCase0():
    p0 = box(0, 0, 1, 1)
    p1 = box(0.5, 0.5, 1.5, 1.5)
    return p0, p1, 'Test case 0'

def TestCase1():
    p0 = box(0, 0, 1, 1)
    p1 = box(1.1, 1.1, 2.1, 2.1)
    return p0, p1, 'Test case 1'

def GetVertices(polygon):
    ring = get_exterior_ring(polygon)
    points = [get_point(ring, i) for i in range(get_num_points(ring))]
    x = [get_x(p) for p in points]
    y = [get_y(p) for p in points]
    return x, y

def MergePolygons(p0, p1):
    return p0

def RunTestCase(testCase):

    p0, p1, title = testCase()

    pm = p0

    x0, y0 = GetVertices(p0)
    x1, y1 = GetVertices(p1)
    xm, ym = GetVertices(pm)

    xx = x0 + x1 + xm
    yy = y0 + y1 + ym
    dx = 0.1*(max(xx) - min(xx))
    dy = 0.1*(max(yy) - min(yy))
    _axis = [min(xx) - dx, max(xx) + dx, min(yy) - dy, max(yy) + dy]

    figure()

    subplot(1, 2, 1)
    plot(x0, y0, '-o')
    plot(x1, y1, '-o')
    axis(_axis)
    subplot(1, 2, 2)
    plot(xm, ym, '--o')
    axis(_axis)
    suptitle(title)

if __name__ == '__main__':

    RunTestCase(TestCase0)
    RunTestCase(TestCase1)

    show()
