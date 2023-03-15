These files are used to generated the corresponding C++ header files.
To regenerate the code, call the FEniCS Form Compiler (ffc):

    ffc -l dolfin Laplacian.ufl
    ffc -l dolfin LinearSpace2D.ufl
