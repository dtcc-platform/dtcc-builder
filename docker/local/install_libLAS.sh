#!/bin/bash
git clone git://github.com/libLAS/libLAS.git && cd libLAS && mkdir build && cd build && cmake .. -DWITH_LASZIP=TRUE && make all -j 4 
sudo cp ./bin/Release/liblas.so* /usr/lib/ && sudo cp ../include/* /usr/local/include/ -r

