#!/usr/bin/env bash
#
# Demo of the mesh generation pipeline, including
#
# 1. Height map generation:  vc-generate-heightmap
# 2. City model generation:  vc-generate-citymodel
# 3. Volume mesh generation: vc-generate-mesh

# Generate height map
../vc-generate-heightmap/build/vc-generate-heightmap Parameters.json

# Generate city model
../vc-generate-citymodel/build/vc-generate-citymodel Parameters.json

# Generate mesh
../vc-generate-mesh/build/vc-generate-mesh Parameters.json
