# Usage

## Overview

DTCC Builder consists of a C++ library and a number of programs
(binaries). The two main programs are

    dtcc-generate-citymodel
    dtcc-generate-mesh

The first of these programs is used to [generate city models from raw
data](#generating-city-models) and the second program is used to
[generate meshes for a city model](#generating-meshes). Both programs
are described in detail below.

A good starting point is to first run the main demo of DTCC Builder by
entering the `demo` directory and issuing the command

    ./dtcc-builder-demo

This will present a menu from which you may select a demo dataset,
after which a city model and a number of output meshes will be created.

The output data may be found in the corresponding subdirectory of the
`data` directory and consist of several data files in JSON and
[Paraview](https://www.paraview.org/) format. Both the data formats
and how to visualize the generated city models and meshes are
described in detail below.

> **Note:** Running `dtcc-generate-mesh` on an Apple M1 processor is 
> currently not fully supported due to problems with the M1 processor and 
> the triangle library. 
 
> **Note:** To run the demo, you first need to build and install DTCC Builder. You
must also download the demo datasets as described in the
[installation instructions](./installation.md) for DTCC Builder.

> **Note:** The demo simply runs the two programs
`dtcc-generate-citymodel` and `dtcc-generate-mesh` in sequence for the
selected dataset.

## Generating city models (`dtcc-generate-citymodel`)

The program `dtcc-generate-citymodel` is used to generate a city model
from a set of point clouds and cadastral data.

### Input data

The following input data are needed:

* **Point cloud data** in LAS/LAZ format consisting of one or more files
  with suffix `.las` or `.laz`.
* **Cadastral data** in [shapefile format](https://en.wikipedia.org/wiki/Shapefile)
  named `PropertyMap.[shp,shx,dbf,prj,cpg]`.
* **Parameters** used to control the city model generation stored
  as a JSON file named `Parameters.json` (optional).

If no command-line argument is given, it is assumed that the current
working directory contains the input data:

    dtcc-generate-citymodel

If a directory is given as command-line argumennt, the given directory
is searched for the input data:

    dtcc-generate-citymodel <path to data directory>

If a parameter file is given as argument, the specified
`DataDirectory` parameter is searched for the input data:

    dtcc-generate-citymodel <path to parameter file>

### Output data

* `CityModel.json` - city model in DTCC JSON format
* `DSM.json` - digital surface map in DTCC JSON format
* `DSM.vts` - digital surface map in VTK structured grid format
* `DTM.json` - digital terrain map in DTCC JSON format
* `DTM.vts` - digital terrain map in VTK structured grid format

In addition, timings and parameters are stored as
`dtcc-generate-citymodel-timings.json` and
`dtcc-generate-citymodel-parameters.json`.

## Generating meshes (`dtcc-generate-mesh`)

The program `dtcc-generate-mesh` is used to generate meshes from a
city model and a digital terrain map.

### Input data

The following input data are needed:

* **City model** in DTCC JSON format named `CityModel.json`.
* **Digital terrain map** in DTCC JSON format named `DTM.json`.
* **Parameters** used to control the mesh generation stored
  as a JSON file named `Parameters.json` (optional).

If no command-line argument is given, it is assumed that the current
working directory contains the input data:

    dtcc-generate-mesh

If a directory is given as command-line argumennt, the given directory
is searched for the input data:

    dtcc-generate-mesh <path to data directory>

If a parameter file is given as argument, the specified
`DataDirectory` parameter is searched for the input data:

    dtcc-generate-mesh <path to parameter file>

### Output data

- `CityModelSimple.json` - simplified city model in DTCC JSON format
- `GroundSurface.json` - surface mesh of ground in DTCC JSON format
- `GroundSurface.vtu` - surface mesh of ground in VTK unstructured grid format
- `BuildingSurface.json` - surface mesh of buildings in DTCC JSON format
- `BuildingSurface.vtu` - surface mesh of buildings in VTK unstructured grid format
- `CitySurface.json` - surface mesh of ground and buildings in DTCC JSON format
- `CitySurface.vtu` - surface mesh of ground and buildings in VTK unstructured grid format
- `CityMesh.json` - volume mesh of city in DTCC JSON format
- `CityMesh.vtu` - volume mesh of city in VTK unstructured grid format

In addition, timings and parameters are stored as
`dtcc-generate-mesh-timings.json` and
`dtcc-generate-mesh-parameters.json`.

## Visualizing results

Generated data files in DTCC JSON format may be opened and visualized
using [DTCC Viewer](https://viewer.dtcc.chalmers.se).

Generated data files in VTK structured/unstructured grid format may be
opened and visualized using [Paraview](https://www.paraview.org/).

## Parameters

DTCC Builder may be controlled using a set of parameters specified in
JSON format. The parameters file may either be supplied as a
command-line argument or stored in a file named `Parameters.json` in
the data directory.

All data files are assumed to be located in a directory determined by
the parameter `DataDirectory`:

    DataDirectory = directory for input data files

Generated data files will be stored in a directory determined by the
parameter `OutputDirectory`:

    OutputDirectory = directory for generated data files

When parsing data from original data files (LAS point clouds and SHP
files), a nonzero origin may be specified to offset the coordinate
system relative to the origin. This has the advantage that very large
values for the coordinates may be avoided (which is good for numerical
stability):

    X0 = x-coordinate of new origin
    Y0 = y-coordinate of new origin

The offset `(X0, Y0)` is subtracted from the original coordinates
during processing. In the simplest case, the offset should be set to
the coordinates of the lower left (south-west) corner of the domain
covered by the data.

Height maps, city models, and meshes are generated for a rectangular
domain with coordinates relative to the new origin specified by `X0`
and `Y0`:

    XMin = x-coordinate for lower left corner
    YMin = y-coordinate for lower left corner
    XMax = x-coordinate for upper right corner
    YMax = y-coordinate for upper right corner

In the simplest case, the lower left corner should be set to `(XMin,
YMin) = (0, 0)` and the upper right corner should be set to `(XMax,
YMax) = (Width, Height)`.

Alternatively, the domain may be determined by the bounding box of the
point cloud(s) by. If `AutoDomain` is `true`, then `XMin`, `YMin`,
`XMax`, `YMax` are automatically determined (and their parameter
values ignored):

    AutoDomain = true/false

When generating elevation models from LAS point cloud data, the
`ElevationModelResolution` parameter determines the resolution of the grid
onto which the height map is sampled:

    ElevationModelResolution = resolution of elevation models

When generating city models from SHP file data, the
`MinimalBuildingDistance` parameter determines a minimal distance
between buildings. Buildings that are closer than the specified
distance are automatically merged to avoid overlapping buildings or
buildings that are very close (which may otherwise upset the mesh
generation):

    MinBuildingDistance = minimal distance between buildings

When generating the volume mesh, the `DomainHeight` parameter
determines the height of the domain relative to the mean ground level:

    DomainHeight = height of computational domain (volume mesh)

When generating both volume and visualization meshes, the
`MeshResolution` parameter determines the maximum size (diameter) of
the mesh cells:

    MeshResolution = resolution of computational mesh (mesh size)

Both volume and visualization meshes may be generated with or without
displacing the ground level outside of buildings. If the `FlatGround`
parameter is set to `true`, then the ground is kept flat:

    FlatGround = true / false

The surface mesh generation produces an additional smoothed version of
the ground surface. The number of smoothing iterations is controlled
by the `GroundSmoothing` parameter:

    GroundSmoothing = number of smoothing iterations

> **Note**: The list of parameters above is only partly complete since
experimental parameters may be added/removed during development. For
a complete list of  parameters, refer to the parameter files
`dtcc-generate-[citymodel,mesh].json` generated by running the demo.


