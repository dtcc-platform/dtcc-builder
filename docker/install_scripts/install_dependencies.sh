#!/bin/bash

apt-get update && apt-get install -y \
    locales \
    sudo \
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
    uuid-dev

./install_assimp.sh
./install_libgeos.sh
./install_libLAS.sh
./install_VTK.sh
