#!/bin/bash

apt-get update && apt-get install -y \
    build-essential \
    cmake\
    nlohmann-json3-dev \
    libshp-dev \
    libpugixml-dev \
    libproj-dev \
    libtriangle-dev \
    libnetcdf-c++4-dev \
    libpng++-dev \
    clang-format \
    clang-tidy \
    doxygen graphviz\
    libgeotiff-dev\
    zlib1g \
    zlib1g-dev \
    zlibc \
    automake\
    colordiff\
    libuuid1 \
    fenics \
    gdal-bin \
    git \
    libgdal-dev \
    nano \
    rsync \
    wget \
    moreutils \
    unzip \
    jq \
    uuid-dev \
    dialog \
    python3 \
    python3-pip \
    python3-dev \
    libprotobuf-dev \
    protobuf-compiler

./install_assimp.sh
./install_libgeos.sh
./install_libLAS.sh
./install_VTK.sh
./install_py_libs.sh