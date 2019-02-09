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

vc-generate-mesh depends on FEniCS version 2018.1 and Triangle.
Use a standard Ubuntu/FEniCS image and run the commands

```
sudo apt-get update
sudo apt-get install libtriangle-dev nlohmann-json-dev
```

### Installing

In preparation.

```
Give examples
```

## Deployment

In preparation.

## Code organization

The code in this repository is organized as a collection of independent but interoperable components. Each component may be implemented using different libraries, and languages (C++, Python, ...) but follows a common naming scheme and provides a standardized command-line interface.

The following list summarizes the currently implemented (and planned components):

* vc-generate-heightmap (PNG to JSON conversion of height maps)
* vc-generate-mesh      (finite element mesh generation)
* vc-simulate-foo       (simulator in preparation)
* vc-simulate-bar       (simulator in preparation)

Common C++ code that is used across components is header only and is placed in the common directory `include`. The common code should have no (or minimal) external dependencies.

## Coding style

VCCore uses Microsoft C# coding style:

```
ClassName
MemberFunction
PublicMemberVariable
PrivateMemberVariable
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
