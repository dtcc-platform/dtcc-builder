#!/usr/bin/bash

python3 -m pip install \
    numpy \
    laspy[lazrs] \
    protobuf==3.20.* \
	pybind11==2.10.* \
	shapely \
    setuptools

apt-get install -y python3-fiona