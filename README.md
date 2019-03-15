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

RH 2000
SWEREF99 12 00 (EPSG:3007)

https://zeus.slu.se
Lantmäteriet:Fastighetskartan bebyggelse vektor EPSG:3006
Lantmäteriet:Laserdata vektor EPSG:3006
by_get

Geodata Extraction Tool” eller GET,

### OpenStreetMap (EPSG:4326)

OpenStreetMap uses the [WGS-84](https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84) coordinate system which is also known as EPSG:4326.

### Lantmäteriet (EPSG:3006)

Lantmäteriet uses the SWEREF99 TM coordinate system which is also known as EPSG:3006.

### VCCore

FIXME: Need to decide on coordinate system. So far we have used
SWEREF99 12 00 (EPSG:3007) which is more accurate for Gothenburg.

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
