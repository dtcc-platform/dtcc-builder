!#/bin/bash
sudo add-apt-repository ppa:fenics-packages/fenics
sudo add-apt-repository ppa:ubuntugis/ppa

sudo apt-get update && apt-get install -y \
 locales \
    sudo \
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
    automake\
    colordiff\
    libuuid1 \
    lftp \
    fenics \
    gdal-bin \
    libgdal-dev \
    git \
    git-lfs \
    dos2unix \
    nano \
    rsync \
    wget \
    moreutils \
    unzip \
    jq

git clone -n https://github.com/LASzip/LASzip.git &  cd LASzip && git checkout 585a940c8d80f039fd1294ecd1411440938d7241 && ./autogen.sh && ./configure && make all -j 4 && make install && mkdir /usr/local/include/laszip/ && mv /usr/local/include/las*.hpp /usr/local/include/laszip/

git clone git://github.com/libLAS/libLAS.git && cd libLAS && mkdir build && cd build && cmake .. -DWITH_LASZIP=TRUE && make all -j 4 && cp ./bin/Release/liblas.so* /usr/lib/ && cp ../include/* /usr/local/include/ -r

