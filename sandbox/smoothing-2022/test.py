# Simple test to compare simple (ad hoc) mesh smoothing with Laplacian smoothing
# Anders Logg 2022-11-30

from fenics import *

# Parameters
n = 32
A = 0.65
a = 0.15
x_0 = 0.75


def smooth_simple(mesh, f):
    # Create new mesh
    new_mesh = Mesh()
    editor = MeshEditor()
    editor.open(new_mesh, "triangle", 2, 2)
    editor.init_vertices(mesh.num_vertices())
    editor.init_cells(mesh.num_cells())
    for i, c in enumerate(mesh.cells()):
        editor.add_cell(i, c)

    # Set new coordinates
    for i, v in enumerate(mesh.coordinates()):
        x, y = v
        y += (1.0 - y) * f(x)
        editor.add_vertex(i, (x, y))

    editor.close()

    return new_mesh


def smooth_laplacian(mesh, f):
    # Boundary value
    class BoundaryValue(UserExpression):
        def eval(self, values, x):
            if near(x[1], 0.0):
                values[0] = f(x[0])
            else:
                values[0] = 0.0

        def value_shape(self):
            return ()

    # Boundary
    boundary = "on_boundary && (near(x[1], 0.0) || near(x[1], 1.0))"

    # Solve PDE for displacement
    parameters["reorder_dofs_serial"] = False
    V = FunctionSpace(mesh, "Lagrange", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = dot(grad(u), grad(v)) * dx
    L = Constant(0) * v * dx
    u = Function(V)
    bv = BoundaryValue()
    bc = DirichletBC(V, bv, boundary)
    solve(a == L, u, bc)

    # Create new mesh
    new_mesh = Mesh()
    editor = MeshEditor()
    editor.open(new_mesh, "triangle", 2, 2)
    editor.init_vertices(mesh.num_vertices())
    editor.init_cells(mesh.num_cells())
    for i, c in enumerate(mesh.cells()):
        editor.add_cell(i, c)

    # Set new coordinates
    dy = u.vector()[:]
    for i, v in enumerate(mesh.coordinates()):
        x, y = v
        y += dy[i]
        editor.add_vertex(i, (x, y))

    editor.close()

    return new_mesh


def smooth_elastic(mesh, f):
    # Boundary value for x-coordinate
    class BoundaryValueX(UserExpression):
        def eval(self, values, x):
            values[0] = 0.0

        def value_shape(self):
            return ()

    # Boundary value for y-coordinate
    class BoundaryValueY(UserExpression):
        def eval(self, values, x):
            if near(x[1], 0.0):
                values[0] = f(x[0])
            else:
                values[0] = 0.0

        def value_shape(self):
            return ()

    # Boundaries
    boundary_x = (
        "on_boundary && (near(x[0], 0.0) || near(x[0], 2.0) || near(x[1], 0.0))"
    )
    boundary_y = "on_boundary && (near(x[1], 0.0) || near(x[1], 1.0))"

    # Define elastic model
    E = 1.0
    nu = 0.499
    lmbda = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))

    def epsilon(u):
        return 0.5 * (grad(u) + grad(u).T)

    def sigma(u):
        return lmbda * div(u) * Identity(2) + 2 * mu * epsilon(u)

    # Solve PDE for displacement
    parameters["reorder_dofs_serial"] = False
    V = VectorFunctionSpace(mesh, "P", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(sigma(u), epsilon(v)) * dx
    L = dot(Constant((0, 0)), v) * dx
    u = Function(V)
    bvx = BoundaryValueX()
    bvy = BoundaryValueY()
    bcx = DirichletBC(V.sub(0), bvx, boundary_x)
    bcy = DirichletBC(V.sub(1), bvy, boundary_y)
    solve(a == L, u, [bcx, bcy])

    # Create new mesh
    new_mesh = Mesh()
    editor = MeshEditor()
    editor.open(new_mesh, "triangle", 2, 2)
    editor.init_vertices(mesh.num_vertices())
    editor.init_cells(mesh.num_cells())
    for i, c in enumerate(mesh.cells()):
        editor.add_cell(i, c)

    # Set new coordinates
    dxy = u.vector()[:]
    N = len(dxy)
    for i, v in enumerate(mesh.coordinates()):
        x, y = v
        x += dxy[i]
        y += dxy[N // 2 + i]
        editor.add_vertex(i, (x, y))

    editor.close()

    return new_mesh


def smooth_hyperelastic(mesh, f):
    # Boundary value for x-coordinate
    class BoundaryValueX(UserExpression):
        def eval(self, values, x):
            values[0] = 0.0

        def value_shape(self):
            return ()

    # Boundary value for y-coordinate
    class BoundaryValueY(UserExpression):
        def eval(self, values, x):
            if near(x[1], 0.0):
                values[0] = f(x[0])
            else:
                values[0] = 0.0

        def value_shape(self):
            return ()

    # Boundaries
    boundary_x = (
        "on_boundary && (near(x[0], 0.0) || near(x[0], 2.0) || near(x[1], 0.0))"
    )
    boundary_y = "on_boundary && (near(x[1], 0.0) || near(x[1], 1.0))"

    # Define elastic model
    E = 1.0
    nu = 0.45  # does not converge with 0.49
    lmbda = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))

    # Solve PDE for displacement
    parameters["reorder_dofs_serial"] = False
    V = VectorFunctionSpace(mesh, "P", 1)
    du = TrialFunction(V)
    v = TestFunction(V)
    u = Function(V)
    I = Identity(2)
    F = I + grad(u)
    C = F.T * F
    Ic = tr(C)
    J = det(F)
    psi = (mu / 2) * (Ic - 3) - mu * ln(J) + (lmbda / 2) * (ln(J)) ** 2
    Pi = psi * dx
    F = derivative(Pi, u, v)
    J = derivative(F, u, du)
    bvx = BoundaryValueX()
    bvy = BoundaryValueY()
    bcx = DirichletBC(V.sub(0), bvx, boundary_x)
    bcy = DirichletBC(V.sub(1), bvy, boundary_y)
    solve(F == 0, u, [bcx, bcy], J=J)

    # Create new mesh
    new_mesh = Mesh()
    editor = MeshEditor()
    editor.open(new_mesh, "triangle", 2, 2)
    editor.init_vertices(mesh.num_vertices())
    editor.init_cells(mesh.num_cells())
    for i, c in enumerate(mesh.cells()):
        editor.add_cell(i, c)

    # Set new coordinates
    dxy = u.vector()[:]
    N = len(dxy)
    for i, v in enumerate(mesh.coordinates()):
        x, y = v
        x += dxy[i]
        y += dxy[N // 2 + i]
        editor.add_vertex(i, (x, y))

    editor.close()

    return new_mesh


def smooth_winslow(mesh, f):
    # Wrapper for boundary value
    class BoundaryValue(UserExpression):
        def eval(self, values, x):
            if near(x[1], 0.0):
                values[0] = 0.0
                values[1] = f(x[0])
            else:
                values[0] = 0.0
                values[1] = 0.0

        def value_shape(self):
            return (2,)

    # Solve PDE for displacement
    parameters["reorder_dofs_serial"] = False
    V = VectorFunctionSpace(mesh, "P", 1)
    u = Function(V)
    v = TestFunction(V)
    alpha = (u[0].dx(1)) ** 2 + (u[1].dx(1)) ** 2
    beta = u[0].dx(0) * u[0].dx(1) + u[1].dx(0) * u[1].dx(1)
    gamma = (u[0].dx(0)) ** 2 + (u[1].dx(0)) ** 2
    A = as_tensor(((alpha, -beta), (-beta, gamma)))
    F = (
        inner(grad(v[0]), dot(A, grad(u[0]))) * dx
        + inner(grad(v[1]), dot(A, grad(u[1]))) * dx
    )
    bv = BoundaryValue()
    bc = DirichletBC(V, bv, "on_boundary && (near(x[1], 0.0) || near(x[1], 1.0))")
    solve(F == 0, u, bc)

    # Create new mesh
    new_mesh = Mesh()
    editor = MeshEditor()
    editor.open(new_mesh, "triangle", 2, 2)
    editor.init_vertices(mesh.num_vertices())
    editor.init_cells(mesh.num_cells())
    for i, c in enumerate(mesh.cells()):
        editor.add_cell(i, c)

    # Set new coordinates
    dxy = u.vector()[:]
    N = len(dxy)
    for i, v in enumerate(mesh.coordinates()):
        x, y = v
        x += dxy[i]
        y += dxy[N // 2 + i]
        editor.add_vertex(i, (x, y))

    editor.close()

    return new_mesh


# Define ground height
f = lambda x: A * exp(-((x - x_0) ** 2) / (2 * a**2))
# f = lambda x : 0.5

# Create mesh
mesh = RectangleMesh(Point(0, 0), Point(2, 1), 2 * n, n)

# Smooth mesh
mesh_simple = smooth_simple(mesh, f)
mesh_laplacian = smooth_laplacian(mesh, f)
mesh_elastic = smooth_elastic(mesh, f)
# mesh_hyperelastic = smooth_hyperelastic(mesh, f)
# mesh_winslow = smooth_winslow(mesh, f)

# Write to file
File("mesh_simple.pvd") << mesh_simple
File("mesh_laplacian.pvd") << mesh_laplacian
File("mesh_elastic.pvd") << mesh_elastic
# File('mesh_hyperelastic.pvd') << mesh_hyperelastic
# File('mesh_winslow.pvd') << mesh_winslow
