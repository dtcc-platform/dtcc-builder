#!/bin/bash

apt-get update && apt-get install -y \
    locales \
    sudo \
	build-essential \
    cmake \
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
    libeigen3-dev 

# get latest cmake 
# wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | sudo tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null
# apt-add-repository "deb https://apt.kitware.com/ubuntu/ $(lsb_release -cs) main"
# apt-get update
# apt-get install -y \
#     cmake 