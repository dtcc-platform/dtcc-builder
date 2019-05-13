# VCCore - VirtualCity@Chalmers Core

VirtualCity@Chalmers is a multidisciplinary research project at
Chalmers University of Technology involving researchers from
mathematics, architecture, civil engineering and computer science. The
aim is to develop an open multimodal data, simulation and
visualization platform for interactive planning, design, exploration,
experimentation and optimization of cities.

This repository contains the server side functionality, including
mesh generation, solvers and data processing.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

libshp-dev liblas-dev liblas-c-dev

#### vc-generate-heightmap

Use a clean Ubuntu 18.04 image with the following packages:

```
sudo apt-get update
sudo apt-get install nlohmann-json-dev libmagick++-dev
```

#### vc-generate-citymodel

Use a clean Ubuntu 18.04 image with the following packages:

```
sudo apt-get update
sudo apt-get install nlohmann-json-dev libpugixml-dev libproj-dev
```

#### vc-randomize-citymodel

Use a clean Ubuntu 18.04 image with the following packages:

```
sudo apt-get update
sudo apt-get install nlohmann-json-dev
```

#### vc-generate-mesh

Use a clean FEniCS 2018.1 image (Ubuntu 18.04) with the following packages:

```
sudo apt-get update
sudo apt-get install libtriangle-dev nlohmann-json-dev
```
#### vc-web-api

Requires latest flask, but should only be used on cloud.virtualcity.chalmers.se.

### Installing

In preparation.

```
Give examples
```

## Deployment

In preparation.

## Data sources

VCCore makes use of the following data sources:

* Building footprints from [OpenStreetMap](https://www.openstreetmap.org)
* Height maps from [Lantmäteriet](https://www.lantmateriet.se/sv/Kartor-och-geografisk-information/Hojddata/Laserdata/laserdata-nh/)

## Coordinate systems

VCCore users meters as a unit of length, relative to the SWEREF99 TM (EPSG:3006) coordinate system.

FIXME: Should we use SWEREF99 12 00 (EPSG:3007)?

FIXME: What about RH 2000?

https://zeus.slu.se
Lantmäteriet:Fastighetskartan bebyggelse vektor EPSG:3006
Lantmäteriet:Laserdata vektor EPSG:3006

## Parameters

VCCore uses the following global parameters, controlled via a JSON file Parameters.json.

When parsing data from original data files (LAS point clouds and SHP files), a nonzero origin may be specified to offset the coordinate system relative to the origin. This has the advantage that very large values for the coordinates may be avoided (which is good for numerical stability).

    X0 - x coordinate of new origin
    Y0 - y coordinate of new origin

Height maps, city models and mesh are generated for a square domain with coordinates relative to the new origin specified by `X0` and `Y0`.

    XMin - x coordinate for lower left corner
    YMin = y coordinate for lower left corner
    XMax = x coordinate for upper right corner
    YMax = y coordinate for upper right corner

When generating the height map from LAS point cloud data, the `HeighMapResolution` parameter determines the resolution of the grid on to which the height map is sampled.

    HeightMapResolution - resolution of height map grid

When generating the city model from SHP file data, the `MinimalBuildingDistance` parameter determines a minimal distance between buildings. Buildings that are closer than the specified distance are automatically merged to avoid overlapping buildings or buildings that are very close (which may otherwise upset the mesh generation).

    MinimalBuildingDistance - minimal distance between buildings

When generating the volume mesh, the `DomainHeight` parameter determines the height of the domain relative to the average ground level. The `MeshResolution` parameter determines the maximum size (diameter) of the mesh cells.

    DomainHeight   - height of computational domain (volume mesh)
    MeshResolution - resolution of computational mesh (mesh size)

## Code organization

The code in this repository is organized as a collection of independent but interoperable components. Each component may be implemented using different libraries, and languages (C++, Python, ...) but follows a common naming scheme and provides a standardized command-line interface.

The following list summarizes the currently implemented (and planned components):

* vc-generate-heightmap  (generate JSON height map data)
* vc-generate-citymodel  (generate JSON city model data from OSM data)
* vc-generate-mesh       (generate FEM mesh from city model and height map)
* vc-generate-mesh-batch (generate FEM meshes from a batch of city models)
* vc-randomize-citymodel (randmize JSON city model)
* vc-plot-mesh           (plot JSON mesh, for testing/debugging)
* vc-plot-citymodel      (plot JSON city model, for testing/debugging)
* vc-simulate-foo        (simulator in preparation)
* vc-simulate-bar        (simulator in preparation)

Common C++ code that is used across components is header only and is placed in the common directory `include`. The common code should have no (or minimal) external dependencies.

## Coding style

VCCore uses Microsoft C# coding style (for both C++ and Python code):

```
ClassName
MemberFunction
PublicMemberVariable
privateMemberVariable
argumentName
variableName
```

## Versioning

VCCore uses [CalVer](https://calver.org/) for versioning.

## Authors (in alphabetical order)

* [Anders Logg](http://anders.logg.org)
* [Vasilis Naserentin](https://www.chalmers.se/en/Staff/Pages/vasnas.aspx)

## License

VCCore is licensed under TBD.

## Acknowledgments
