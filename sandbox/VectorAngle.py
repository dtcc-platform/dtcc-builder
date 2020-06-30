# Prototyping function for cheap computation of "angle" between vectors
# Anders Logg 2020-06-24

from pylab import *

class Vector:
    def __init__(self, vx, vy):
        self.x = vx
        self.y = vy

def VectorAngle2D(u, v):
    "Return strictly increasing function of angle of v relative to u"
    u2 = u.x*u.x + u.y*u.y
    v2 = v.x*v.x + v.y*v.y
    sin = u.x*v.y - u.y*v.x
    cos = u.x*v.x + u.y*v.y
    a = sin*sin / (u2*v2)
    if sin > 0.0:
        return a if cos > 0.0 else 2.0 - a
    else:
        return -a if cos > 0.0 else a - 2.0

a0 = 3.4
ur = 3.0
vr = 2.5
x = linspace(-pi, pi, 1000)
y = []
for a in x:
    u = Vector(ur*cos(a0), ur*sin(a0))
    v = Vector(vr*cos(a0 + a), vr*sin(a0 + a))
    y.append(VectorAngle2D(u, v))

plot(x, y)
show()
