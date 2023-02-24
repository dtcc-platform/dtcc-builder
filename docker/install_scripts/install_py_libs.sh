#!/usr/bin/bash

python3 -m pip install \
    numpy>=1.23.* \
    laspy[lazrs] \
    protobuf==3.20.* \
	pybind11==2.10.* \
	shapely \
    setuptools \
    meshio

apt-get install -y python3-fiona python3-rasterio