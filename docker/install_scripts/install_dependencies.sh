#!/bin/bash
#add-apt-repository -y ppa:ubuntugis/ubuntugis-unstable
#apt-get update

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
    libpthread-stubs0-dev \
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
    libassimp-dev \
    protobuf-compiler

#these seem necessary for the python bindings to work
if [ -f /usr/lib/x86_64-linux-gnu/libassimp.so.5 ]; then
    ln -s /usr/lib/x86_64-linux-gnu/libassimp.so.5 /usr/lib/
fi
if [ -f /usr/lib/aarch64-linux-gnu/libassimp.so.5]; then
    ln -s /usr/lib/aarch64-linux-gnu/libassimp.so.5 /usr/lib/
fi

# remove old version of libgeos and install new version from source
apt remove -y libgeos-dev
./install_libgeos.sh

# ./install_assimp.sh
# ./install_libLAS.sh
./install_VTK.sh
#./install_py_libs.sh
#./install_dtccmodel.sh
