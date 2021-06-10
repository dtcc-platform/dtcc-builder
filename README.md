# Digital Twin Cities Platform - Core

_This README was last revised on 3rd of March 2021_

The Digital Twin Cities Platform is an open-source platform for the
exploration of digital twins for cities. The platform is developed and
maintained by the Digital Twin Cities Centre (DTCC) hosted by Chalmers
University of Technology. The aim is to develop an open multimodal
data, modeling, simulation and visualization platform for interactive
planning, design, exploration, experimentation, and optimization of
cities.

This repository (Core) provides software for data processing,
modeling, and simulation.

![](images/hammarkullen.jpg)

*The image shows a mesh of the Hammarkullen area in Gothenburg
 generated from public data from the [Swedish mapping, cadastral and
 land registration authority (Lantmäteriet)](https://www.lantmateriet.se/).*

## Getting started

These instructions will get you up and running on your local machine
for development and testing purposes.

### Downloading the software

To download the software, clone the repository by the following command:

    git clone https://gitlab.com/dtcc3d/core.git

Alternatively, you may want to use the SSH protocol:

    git clone git@gitlab.com:dtcc3d/core.git

This will create a directory named `core` containing the full source code.

**Note:** If you are using Windows, you might first want to make sure
that Git does not convert Unix-style file endings on checkout. This
can be accomplished by:

    git config --global core.autocrlf false

### Building the Docker container

The most convenient way to work with the Digital Twin Cities Platform
is to use the custom [Docker](https://www.docker.com/) image, which
contains all the dependencies needed for developing, building, and
running the software.

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
container (virtual machine) in which to run the Digital Twin Cities
Platform:

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

To download demo data, there is an open demo dataset (public) and one available only for development/research purposes. Enter the `data` directory and issue the following command:

    ./dtcc-download-demo-data-public

For the typical new user, this should be done inside the container. Obtaining the public dataset requires no credentials.

For accesing the private demo datasets you will need a set of credentials supplied by the development team. Contact us for more info.

### Building and installing

To build core, use a standard out-of-source CMake build by issuing the following commands from the top level directory:

    mkdir build
    cd build
    cmake ..
    make
    make install

This will build and install all programs into the top level `bin` directory.

### Running the demo

To run a simple demo, enter the `demo` directory and issue the following command:

    ./dtcc-demo

This will generate a set of meshes for the demo data downloaded. The input data as well as the output data generated are stored in the `data` directory. To change dataset or parameters for the demo, edit the script `dtcc-demo` to select a different paremeter set. See below for a description of the parameters.

## Data

### Data sources

The Digital Twin Cities Platform makes use of the following data sources:

* Point clouds (Lantmäteriet:Laserdata vektor; EPSG:3006)
* Property maps (Lantmäteriet:Fastighetskartan bebyggelse vektor; EPSG:3006)
* Road network maps (Lantmäteriet:Vägkartan; EPSG:3006)

Chalmers has a license for downloading data from `http://zeus.slu.se`.

Point cloud data comes in the form of a number square grids large enough to cover the requested domain. Each point cloud is compressed as a RAR file. Uncompress it to get the LAS point cloud file, for example:

    unrar e 09B008_64050_3225_25.rar

This will create the file 09B008_64050_3225_25.las. The lower left corner will in this example be at EPSG:3006 coordinates (6405000, 322500).

Property map data and the road network comes in the form of SHP files (with corresponding SHX, DBF and PRJ files). The files of interest are the ones named `by_get.*` and `vl_*` respectivelly.

### Data formats

WIP: Describe DTCC JSON format.

WIP: Describe CityJSON format.

### Coordinate system

Core uses meters as a unit of length, relative to the SWEREF99 TM (EPSG:3006) coordinate system.

## Parameters

The Digital Twin Cities Platform uses the following global parameters, controlled via a JSON file `Parameters.json`.

All data files are assumed to be located in a directory determined by the
parameter `DataDirectory`. Any generated data files will be stored in the
same location.

    DataDirectory = directory for input/output

When parsing data from original data files (LAS point clouds and SHP files), a nonzero origin may be specified to offset the coordinate system relative to the origin. This has the advantage that very large values for the coordinates may be avoided (which is good for numerical stability).

    X0 = x-coordinate of new origin
    Y0 = y-coordinate of new origin

In other words, the offset `(X0, Y0)` is subtracted from the original coordinates during processing. In the simplest case, the offset should be set to the coordinates of the lower left (south-east) corner of the domain covered by the data.

Height maps, city models and meshes are generated for a rectangular domain with coordinates relative to the new origin specified by `X0` and `Y0`.

    XMin = x-coordinate for lower left corner
    YMin = y-coordinate for lower left corner
    XMax = x-coordinate for upper right corner
    YMax = y-coordinate for upper right corner

In the simplest case, the lower left corner should be set to `(XMin, YMin) = (0, 0)` and the upper right corner should be set to `(XMax, YMax) = (Width, Height)`.

Alternatively, the domain may be determined by the bounding box of the point cloud(s) by. If `AutoDomain` is `true`, then `XMin`, `YMin`, `XMax`, `YMax` are automatically determined (and their parameter values ignored).

    AutoDomain = true/false

When generating the height map from LAS point cloud data, the `HeighMapResolution` parameter determines the resolution of the grid onto which the height map is sampled.

    HeightMapResolution = resolution of height map grid

When generating the city model from SHP file data, the `MinimalBuildingDistance` parameter determines a minimal distance between buildings. Buildings that are closer than the specified distance are automatically merged to avoid overlapping buildings or buildings that are very close (which may otherwise upset the mesh generation).

    MinimalBuildingDistance = minimal distance between buildings

When generating the volume mesh, the `DomainHeight` parameter determines the height of the domain relative to the average ground level.

    DomainHeight = height of computational domain (volume mesh)

When generating both volume and visualization meshes, the `MeshResolution` parameter determines the maximum size (diameter) of the mesh cells.

    MeshResolution = resolution of computational mesh (mesh size)

Both volume and visualization meshes may be generated with our without displacing the ground level outside of buildings. If the `FlatGround` parameter is set to `true`, then the ground is kept flat.

    FlatGround = true / false

The surface mesh generation produces an additional smoothed version of the ground surface. The number of smoothing iterations is controlled by the `GroundSmoothing` parameter.

    GroundSmoothing = number of smoothing iterations

## Design

### Code organization

The Digital Twin Cities Platform is organized as a collection of independent but interoperable components. Each component may be implemented using different libraries, and languages (C++, Python, ...) but follows a common naming scheme and provides a standardized command-line interface.

Common C++ code that is used across components is *header only* and is placed in the common directory `include`. The common code should have no (or minimal) external dependencies.

### Coding style

The Digital Twin Cities Platform uses Microsoft C# coding style (for both C++ and Python code):

```
ClassName
MemberFunction
PublicMemberVariable
privateMemberVariable
argumentName
variableName
```

Code formatting is enforced using [ClangFormat](https://clang.llvm.org/docs/ClangFormat.html) as defined in the top-level `.clang-format` file. The style is based on the default LLVM style with minimal modifications.

Algorithms should be implemented as static functions in separate classes (for example in the class `MeshGenerator` rather than in the class `Mesh`). This means that pure data classes (like `Mesh`) can be kept clean with only class data and functions for data access and initialization.

### Versioning

The Digital Twin Cities Platform uses [CalVer](https://calver.org/) for versioning.

## Authors (in order of appearance)

* [Anders Logg](http://anders.logg.org)
* [Vasilis Naserentin](https://www.chalmers.se/en/Staff/Pages/vasnas.aspx)
* [Dag Wästerberg](http://www.ramboll.se)
* [Orfeas Eleutheriou](http://orfeasel.com/)
* [Anton Olsson](mailto:anton.j.olsson@bredband.net)
* [Anton Annlöv](mailto:annlova@student.chalmers.se)

Part of this code is contributed by ReSpace AB under the MIT License.

## Output

- DSM.json - DSM generated from point cloud (GridField2D)
- DTM.json - DTM generated from point cloud (GridField2D)
- DSM.vts - same as above in Paraview .vts (structured VTK)
- DTM.vts - same as above in Paraview .vts (structured VTK)
- CityModelRaw.json    - city model based on unprocessed shapefile data
- CityModelClean.json  - city model based on cleaned data (fix orientation, remove duplicate vertices)
- CityModelSimple.json - city model generated by simplifying CityModelClean (merged buildings)
- BuildingMesh.json - surface mesh of all buildings generated from CityModelClean (Surface3D)
- GroundMesh.json   - surface mesh of ground generated from DTM (Surface3D)
- BuildingMesh.vtu - same as above in Paraview .vtu (unstructured VTK)
- GroundMesh.vtu   - same as above in Paraview .vtu (unstructured VTK)
- Mesh.json     - volume mesh of whole domain generated from CityModelSimple (Mesh3D)
- Boundary.json - boundary of Mesh (Surface3D)
- Step[31-35][Mesh/Boundary].vtu - debugging output, Step35 = final Mesh


## License

The Digital Twin Cities Platform is licensed under the [MIT license](https://opensource.org/licenses/MIT).
Copyright is held by the individual authors as listed at the top of each source file.

## Acknowledgments

This work is part of the Digital Twin Cities Centre supported by Sweden’s Innovation Agency VINNOVA under Grant No. XXX.

WIP: Add grant number.
