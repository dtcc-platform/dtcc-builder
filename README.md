# DTCC Builder

_This README was last revised on 22nd of March 2022_

The DTCC Platform is an open-source platform for the exploration of
digital twins for cities. The platform is developed and maintained by
the Digital Twin Cities Centre (DTCC) hosted by Chalmers University of
Technology. The aim is to develop an open modelling, simulation and
visualisation platform for interactive planning, design, and
exploration of cities.

This repository, DTCC Builder, provides software for generation of 3D
models (meshes) from raw data.

![](images/hammarkullen.jpg)

*The image shows a mesh of the Hammarkullen area in Gothenburg
 generated from public data from the [Swedish mapping, cadastral and
 land registration authority (Lantmäteriet)](https://www.lantmateriet.se/).*

## Getting started

These instructions will get you up and running on your local machine
for development and testing purposes.

### Downloading the software

To download the software, clone the repository by the following command:

    git clone https://gitlab.com/dtcc-platform/dtcc-builder.git -b develop

Alternatively, you may want to use the SSH protocol:

    git clone git@gitlab.com:dtcc-platform/dtcc-builder.git

This will create a directory named `dtcc-builder` containing the full
source code.

**Note:** If you are using Windows, you might first want to make sure
that Git does not convert Unix-style file endings on checkout. This
can be accomplished by:

    git config --global dtcc-builder.autocrlf false

### Building the Docker container

The most convenient way to work with the DTCC Platform is to use the
custom [Docker](https://www.docker.com/) image, which contains all the
dependencies needed for developing, building, and running the
software.

The first step is to download and install
[Docker](https://www.docker.com/). After the installation is complete,
continue with the following steps.

On Linux/MacOS, enter the `docker` directory and issue the following
command:

    ./dtcc-pull-image

On Windows, you should instead enter the `docker/Windows` subdirectory
and issue the following command:

    ./dtcc-pull-image.bat

This creates a Docker image named `dtccimage`.

**Note:** For expert users and/or debugging purposes you can create
the image from scratch by issuing `./dtcc-build-image` and
`dtcc-build-image.bat` on Linux/MacOS and Windows respectively.

Then issue the following commands to create and start a persistent
container (virtual machine) in which to run the DTCC Platform:

    ./dtcc-create-container
    ./dtcc-start-container

On Windows, you should instead issue the following commands:

    ./dtcc-create-container.bat
    ./dtcc-start-container.bat

After completing this step, you should now be inside the Docker
container named `dtcc` and ready to go.

For removing the image and containers, you can use

    ./dtcc-uninstall
    ./dtcc-uninstall.bat

in Linux/MacOS and Windows respectively.

**Note:** The source tree is automatically shared into the Docker
  container. It is recommended that you edit the sources, run Git
  commands, and visualize data *outside* of the Docker container (on
  your native operating system), while building and running the code
  *inside* the Docker container.

### Downloading demo data

To download demo data, there is an open demo dataset (public) and one
available only for development/research purposes. Enter the `data`
directory and issue the following command:

    ./dtcc-download-demo-data-public

For the typical new user, this should be done inside the
container. Obtaining the public dataset requires no credentials.

The public demo data currently contains the following two datasets:

* `Eman2021`: A residential area in the vicinity of the river Emån in
southern Sweden (Småland)
* `Majorna2021`: A residential area within the City of Gothenburg in
  Sweden

The public demo data is kindly provided by Lantmäteriet and published at

    ftp://download-open.lantmateriet.se/GEO-demoomrade

For accesing the private demo datasets you will need a set of
credentials supplied by the development team. Contact us for more
info.

### Building and installing

To build, use a standard out-of-source CMake build by issuing the
following commands from the top level directory:

    mkdir build
    cd build
    cmake ..
    make -j
	make install

This will build and install all programs into the top level `bin`
directory.

**Note:** These commands should be issued *inside* the container.

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

## Design

### Code organization

The DTCC Platform is organized as a collection of independent but
interoperable components. Each component may be implemented using
different libraries, and languages (C++, Python, ...) but follows a
common naming scheme and provides a standardized command-line
interface.

Common C++ code that is used across components is *header only* and is
placed in the common directory `include`. The common code should have
no (or minimal) external dependencies.

### Coding style

The DTCC Platform uses Microsoft C# coding style (for both C++ and
Python code):

```
ClassName
MemberFunction
PublicMemberVariable
privateMemberVariable
argumentName
variableName
```

Code formatting is enforced using
[ClangFormat](https://clang.llvm.org/docs/ClangFormat.html) as defined
in the top-level `.clang-format` file. The style is based on the
default LLVM style with minimal modifications.

Algorithms should be implemented as static functions in separate
classes (for example in the class `MeshGenerator` rather than in the
class `Mesh`). This means that pure data classes (like `Mesh`) can be
kept clean with only class data and functions for data access and
initialization.

### Versioning

The DTCC Platform uses [CalVer](https://calver.org/) for versioning.

## Authors (in order of appearance)

* [Anders Logg](http://anders.logg.org)
* [Vasilis Naserentin](https://www.chalmers.se/en/Staff/Pages/vasnas.aspx)
* [Dag Wästerberg](https://chalmersindustriteknik.se/sv/medarbetare/dag-wastberg/)
* [Orfeas Eleutheriou](http://orfeasel.com/)
* [Anton Olsson](mailto:anton.j.olsson@bredband.net)
* [Anton Annlöv](mailto:annlova@student.chalmers.se)

Part of this code is contributed by ReSpace AB under the MIT License.

## License

The DTCC Platform is licensed under the [MIT
license](https://opensource.org/licenses/MIT).

Copyright is held by the individual authors as listed at the top of
each source file.

## Acknowledgments

This work is part of the Digital Twin Cities Centre supported by
Sweden’s Innovation Agency Vinnova under Grant No. 2019-421 00041.
