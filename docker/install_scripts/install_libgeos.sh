#!/bin/bash

wget http://download.osgeo.org/geos/geos-3.10.2.tar.bz2
tar jxfv geos-3.10.2.tar.bz2
cd geos-3.10.2
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local ..
make -j4
sudo make install
cd ../..
rm -rf geos-3.10.2
rm geos-3.10.2.tar.bz2
