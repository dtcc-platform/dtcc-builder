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
sudo apt-get install nlohmann-json-dev libpugixml-dev
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

## Data

VCCore uses

### Coordinate systems

RT90

WGS-84

SWEREF99 12 00 (EPSG:3007) used in VC
https://epsg.io/3007
more accurate for Gothenburg

SWEREF99 TM (EPSG:3006)



## Code organization

The code in this repository is organized as a collection of independent but interoperable components. Each component may be implemented using different libraries, and languages (C++, Python, ...) but follows a common naming scheme and provides a standardized command-line interface.

The following list summarizes the currently implemented (and planned components):

* vc-generate-heightmap  (generate JSON height map data)
* vc-generate-citymodel  (generate JSON city model data from OSM data)
* vc-generate-mesh       (generate FEM mesh from city model and height map)
* vc-generate-mesh-batch (generate FEM meshes from a batch of city models)
* vc-randomize-citymodel (randmize JSON city model)
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
