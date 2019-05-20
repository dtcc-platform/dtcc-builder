# This file builds the Docker image for VCCore

# Use Phusion base image (minimal Docker-friendly Ubuntu)
FROM phusion/baseimage:0.11

# Install basic dependencies
RUN apt-get update && apt-get install -y \
    nlohmann-json-dev \
    libshp-dev \
    liblas-dev \
    liblas-c-dev \
    libpugixml-dev \
    libproj-dev \
    libtriangle-dev

# Install FEniCS
RUN add-apt-repository ppa:fenics-packages/fenics
RUN apt-get update && apt-get install -y fenics
