## Usage

### Running the demo

To run a simple demo, enter the `demo` directory and issue the
following commands:

    ./demo-generate-citymodel
    ./demo-generate-mesh

This will generate a set of meshes for the demo data downloaded. The
input data as well as the output data generated are stored in the
`data` directory. To change dataset or parameters for the demo, edit
the scripts above to select a different parameter set. See below for a
description of the parameters.

**Note:** These commands should be issued *inside* the container.

### Visualizing results

To visualize the generated meshes, install and open
[Paraview](https://www.paraview.org/). This will allow
you to open and visualize the generated `.vts` and `.vtu`
files in the `data` directory.

## Data

### Data sources

The DTCC Platform makes use of the following data sources:

* Point clouds (Lantmäteriet:Laserdata NH 2019 (laz); EPSG:3006)
* Property maps (Lantmäteriet:Fastighetskartan Bebyggelse; EPSG:3006)
* Road network maps (Lantmäteriet:Vägkartan; EPSG:3006)

Chalmers has a license for downloading data from `http://zeus.slu.se`.

Point cloud data comes in the form of a number square grids large
enough to cover the requested domain. Each point cloud is compressed
as a RAR file. Uncompress it to get the LAS point cloud file, for
example:

    unrar e 09B008_64050_3225_25.rar

This will create the file 09B008_64050_3225_25.las. The lower left
corner will in this example be at EPSG:3006 coordinates (6405000,
322500).

Property map data and the road network comes in the form of SHP files
(with corresponding SHX, DBF and PRJ files). The files of interest are
the ones named `by_get.*` and `vl_*` respectivelly.

### Road network data

Like property maps, road maps also use the SHP format (as well as same
formats for corresponding files). By using the program
`dtcc-generate-roadnetwork`, shapefiles can be converted into JSON
format. From the `bin` directory, the following command is then used:

    ./dtcc-generate-roadnetwork <path to SHP file>

To view the created JSON file, the program `dtcc-plot` can be
used. For usage, see the source file (`dtcc-plot/dtcc-plot`). Plotting
must be done outside the container, using your native system. The
Python package `matplotlib` needs to be installed: `pip3 install
matplotlib` on Mac and Windows.

### Data formats

WIP: Describe DTCC JSON format.

WIP: Describe CityJSON format.

### Coordinate system

The unit of length is metres relative to the SWEREF99 TM (EPSG:3006)
coordinate system.

## Parameters

The DTCC Platform uses the following global parameters,
controlled via a JSON file `Parameters.json`.

All data files are assumed to be located in a directory determined by
the parameter `DataDirectory`. Any generated data files will be stored
in the same location.

    DataDirectory = directory for input/output

When parsing data from original data files (LAS point clouds and SHP
files), a nonzero origin may be specified to offset the coordinate
system relative to the origin. This has the advantage that very large
values for the coordinates may be avoided (which is good for numerical
stability).

    X0 = x-coordinate of new origin
    Y0 = y-coordinate of new origin

In other words, the offset `(X0, Y0)` is subtracted from the original
coordinates during processing. In the simplest case, the offset should
be set to the coordinates of the lower left (south-west) corner of the
domain covered by the data.

Height maps, city models and meshes are generated for a rectangular
domain with coordinates relative to the new origin specified by `X0`
and `Y0`.

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
values ignored).

    AutoDomain = true/false

**Note**: The `AutoDomain` parameter has been temporarily disabled.

When generating elevation models from LAS point cloud data, the
`ElevationModelResolution` parameter determines the resolution of the grid
onto which the height map is sampled.

    ElevationModelResolution = resolution of elevation models

When generating the city model from SHP file data, the
`MinimalBuildingDistance` parameter determines a minimal distance
between buildings. Buildings that are closer than the specified
distance are automatically merged to avoid overlapping buildings or
buildings that are very close (which may otherwise upset the mesh
generation).

    MinBuildingDistance = minimal distance between buildings

When generating the volume mesh, the `DomainHeight` parameter
determines the height of the domain relative to the average ground
level.

    DomainHeight = height of computational domain (volume mesh)

When generating both volume and visualization meshes, the
`MeshResolution` parameter determines the maximum size (diameter) of
the mesh cells.

    MeshResolution = resolution of computational mesh (mesh size)

Both volume and visualization meshes may be generated with or without
displacing the ground level outside of buildings. If the `FlatGround`
parameter is set to `true`, then the ground is kept flat.

    FlatGround = true / false

The surface mesh generation produces an additional smoothed version of
the ground surface. The number of smoothing iterations is controlled
by the `GroundSmoothing` parameter.

    GroundSmoothing = number of smoothing iterations

**Note**: The list of parameters above is only partly complete since
experimental parameters may be added/removed during development. For
the latest list of parameters, refer to the parameter files for the
demos, for example `Majorna2021.json`

### Output

## JSON files

- `DSM.json` - Digital surface map generated from point cloud (`GridField2D`)
- `DTM.json` - Digital terrain map generated from point cloud (`GridField2D`)
- `CityModel.json` - city model generated from property map and point cloud (`CityModel`)
- `CityModelSimple.json` - simplified (merged) city model (`CityModel`)
- `GroundSurface.json` - surface mesh of ground generated from DTM (`Surface3D`)
- `BuildingSurface.json` - surface mesh of all buildings generated from city model (`Surface3D`)
- `CityMesh.json` - volume mesh of city model (`Mesh3D`)
- `CitySurface.json` - surface mesh of city model (`Surface3D`)

## VTS files (structured VTK meshes)

- `DSM.vts` - Digital surface map generated from point cloud
- `DTM.vts` - Digital terrain map generated from point cloud

## VTU files (unstructured VTK meshes)

- `GroundSurface.json` - surface mesh of ground generated from DTM
- `BuildingSurface.json` - surface mesh of all buildings generated from city model
- `CityMesh.json` - volume mesh of city model
- `CitySurface.json` - surface mesh of city model
- `Step[31-35][Mesh/Boundary].vtu` - mesh generation debugging output (intermediate steps with Step35 = final mesh)

**Note**: Some of these data files are only generated when the
parameter `Debug` is set.
