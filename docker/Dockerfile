# This file builds the Docker image for the Digital Twin Cities Platform

# Use Phusion base image (minimal Docker-friendly Ubuntu)
FROM phusion/baseimage:0.11

# Set some variables
ENV USER dtcc
ENV HOME /home/$USER
ENV DIR core

# Avoid stall during installation of tzdata
#ENV TZ=Europe/Stockholm
#RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Update and upgrade
RUN apt-get update && apt-get upgrade -y

# Install basic dependencies
RUN apt-get -y install locales sudo

# Install libraries used by DTCCore
RUN apt-get update && apt-get install -y \
    nlohmann-json-dev \
    libshp-dev \
#    liblas-dev \   #need to build from source to support LASZip
#    liblas-c-dev \
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
    colordiff

# Install FEniCS
RUN add-apt-repository ppa:fenics-packages/fenics
RUN apt-get update && apt-get install -y fenics

# Install GDAL
RUN add-apt-repository ppa:ubuntugis/ppa
RUN apt-get update && apt-get install -y gdal-bin

# Install git with lfs support, dos2unix, nano simple editor, rsync and wget
RUN apt-get update && apt-get install -y git git-lfs dos2unix nano rsync wget

# Get and compile LASZip
RUN git clone -n https://github.com/LASzip/LASzip.git
RUN cd LASzip && git checkout 585a940c8d80f039fd1294ecd1411440938d7241 && ./autogen.sh && ./configure && make all -j 4 && make install
RUN mkdir /usr/local/include/laszip/ && mv /usr/local/include/las*.hpp /usr/local/include/laszip/

# Get and compile libLAS with LASZip suppport

#RUN git clone git://github.com/libLAS/libLAS.git && cd libLAS && mkdir build && cd build && cmake .. -DWITH_LASZIP=TRUE && make all -j 4 && make install
#RUN ln -s /usr/local/lib/liblas.so.3 /usr/lib && ln -s /usr/local/lib/liblaszip.so.4 /usr/lib
#Temp fix since make install has broken since 5/5/2020
RUN git clone git://github.com/libLAS/libLAS.git && cd libLAS && mkdir build && cd build && cmake .. -DWITH_LASZIP=TRUE && make all -j 4 && cp ./bin/Release/liblas.so* /usr/lib/ && cp ../include/* /usr/local/include/ -r

# Get and compile VTK 7.1

RUN wget https://www.vtk.org/files/release/7.1/VTK-7.1.1.tar.gz && tar xvzf VTK-7.1.1.tar.gz && cd VTK-7.1.1/ && mkdir build && cd build && cmake -DVTK_Group_Rendering=OFF -DVTK_BUILD_ALL_MODULES_FOR_TESTS:BOOL=OFF -DVTK_Group_StandAlone=OFF -DModule_vtkCommonCore:BOOL=ON .. && make all -j 4 && make install

# Get and compile assimp latest (27th May 2020)

RUN git clone https://github.com/assimp/assimp && cd assimp && git checkout 1427e67b54906419e9f83cc8625e2207fbb0fcd5 && mkdir build && cd build && cmake .. && make all -j 4 && make install && ln -s /assimp/build/bin/libassimp.so.5 /usr/lib/x86_64-linux-gnu/

# Add user and change to user
RUN useradd -m $USER -G sudo && \
    echo "$USER:$USER" | chpasswd && \
    echo "$USER ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers
USER $USER

# Create shared volume
VOLUME $HOME/$DIR
WORKDIR $HOME/$DIR

# Generate welcome message printed at login
COPY Welcome $HOME/.welcome
RUN echo "cat $HOME/.welcome" >> $HOME/.bashrc

# Start bash login shell
ENTRYPOINT ["/bin/bash", "-l", "-c"]
CMD ["/bin/bash", "-i"]
