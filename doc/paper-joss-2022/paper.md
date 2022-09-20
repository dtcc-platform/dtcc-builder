---
title: 'DTCC Builder: A mesh generator for automatic, efficient, and robust mesh generation for large-scale city modeling and simulation'
tags:
  - mesh generation
  - point cloud
  - cadastral map
  - digital twin
  - city modeling
  - city simulation
  - C++
authors:
  - name: Anders Logg
    orcid: 0000-0002-1547-4773
    equal-contrib: true
    affiliation: "1, 3"
  - name: Vasilis Naserentin
    equal-contrib: true
    affiliation: "1, 3"
  - name: Dag Wästberg
    equal-contrib: true
    affiliation: "2, 3"
affiliations:
 - name: Chalmers University of Technology
   index: 1
 - name: FIXME CIT
   index: 2
 - name: Digital Twin Cities Centre
   index: 3
date: 20 September 2022
bibliography: paper.bib
---

# Summary

DTCC Builder is a mesh generator for automatic, efficient, and robust
mesh generation for large-scale city modeling and simulation.  Using
standard and widely available raw data sources in the form of point
clouds and cadastral data, DTCC Builder generates high-quality 3D
surface and volume meshes, suitable for both visualization and
simulation. In particular, DTCC Builder is capable of generating
large-scale, conforming tetrahedral volume meshes of cities suitable
for finite element (FEM) simulation.

# Statement of need

The interest in creating digital twins, i.e., models of physical
systems that mirror the physical systems in real-time and enable
analysis and prediction, has been rapidly increasing in recent years.
In particular, there has been a surge in the interest for creating
digital twins of cities [@ketzlerDigitalTwinsCities2020]. The creation
of a digital twin of a city involves the creation of a 3D model. Such
3D models may either be created manually, semi-automatically, or in a
fully automatic way from available raw data, often in the form of
point clouds obtained from aerial scanning and cadastral data
(property maps).

3D mesh generation is a very challenging process, especially in the
face of bad quality and low resolution data, which is often the case
for publicly available data for cities. Furthermore, if the 3D meshes
are to be used for modeling and simulation, certain requirements are
posed on the quality of the meshes. DTCC Builder aims to solve these
challenges by automating the mesh generation process in a both robust
and efficient way.

DTCC Builder is part of the open-source digital twin platform
[DTCC Platform](https://platform.dtcc.chalmers.se) developed at the
[Digital Twin Cities Centre](https://dtcc.chalmers.se).

![](demo-majorna.jpg)
*Surface mesh of an area (Majorna) in Gothenburg, generated with DTCC Buider.*

![](demo-majorna-zoom.jpg)
*Detail of surface mesh of an area (Majorna) in Gothenburg, generated with DTCC Builder.*

# Method and implementation

DTCC Builder uses a novel algorithm for mesh generation. The key idea
is to utilize the special geometry of city models to reduce the 3D
mesh generation problem to a 2D problem. A 2D mesh respecting the
polygonal footprints of buildings is generated and then layered to
create 3D mesh. Building heights and ground height are incorporated by
a PDE-based smoothing process. The method and algorihtms are described
in detail in the paper [naserentinXXX2022].

DTCC Builder is implemented in C++ and makes use of several
open-source packages, notably
FEniCS [loggAutomatedSolutionDifferential2012] for solving PDEs,
Triangle [shewchukTriangleEngineering2D1996] for 2D mesh
generation, and
GEOS [geoscontributorsGEOSCoordinateTransformation2021] for
geometric operations.

# Documentation

The documentation for DTCC Builder is published on the
[DTCC GitLab pages](https://gitlab.com/dtcc-platform/dtcc-builder)
as well as on the documentation pages for
[DTCC Platform](https://dtcc.chalmers.se).

# Future work

DTCC Builder currently only provides a C++ and command-line
interface. Future versions will provide a Python interface and also an
online interface as part of
[DTCC Platform](https://platform.dtcc.chalmers.se).

DTCC Builder currently generates city models in Level of Detail (LoD)
1.2 but ongoing work seeks to extend DTCC Bulder to LoD1.3 and LoD2.x.

# Acknowledgements

This work is part of the Digital Twin Cities Centre supported by
Sweden’s Innovation Agency Vinnova under Grant No.  2019-00041.

# References
