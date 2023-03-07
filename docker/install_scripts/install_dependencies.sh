#!/bin/bash


DEBIAN_FRONTEND=noninteractive apt-get update && apt-get install -y \
    build-essential \
    git \
    cmake \
    libpugixml-dev \
    libtriangle-dev \
    nlohmann-json3-dev \
    libpthread-stubs0-dev \
    libuuid1 \
    libeigen3-dev \
    libgdal-dev \
    python3 \
    python3-dolfin \
    python3-pip \
    python3-dev \
    python3-venv \
    libprotobuf-dev \
    libassimp-dev \
    protobuf-compiler \
    libgeos-dev \
    nano \
    rsync \
    wget \
    unzip \
    colordiff \
    clang-format \
    clang-tidy 


# apt-get install -y python3-fiona python3-rasterio

#these seem necessary for the python bindings to work
if [ -f /usr/lib/x86_64-linux-gnu/libassimp.so.5 ]; then
    ln -s /usr/lib/x86_64-linux-gnu/libassimp.so.5 /usr/lib/
fi
if [ -f /usr/lib/aarch64-linux-gnu/libassimp.so.5]; then
    ln -s /usr/lib/aarch64-linux-gnu/libassimp.so.5 /usr/lib/
fi

# remove old version of libgeos and install new version from source
# apt remove -y libgeos-dev
# ./install_libgeos.sh

# ./install_assimp.sh
# ./install_libLAS.sh
#./install_VTK.sh
#./install_py_libs.sh
./install_dtccmodel.sh
