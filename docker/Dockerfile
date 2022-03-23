# This file builds the Docker image for the Digital Twin Cities Platform

# Use Phusion base image (minimal Docker-friendly Ubuntu)
FROM phusion/baseimage:master as base

# Set some variables
ENV USER dtcc
ENV HOME /home/$USER
ENV DIR core

# Avoid stall during installation of tzdata
#ENV TZ=Europe/Stockholm
#RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Add fenics and gdal repos
#RUN add-apt-repository ppa:fenics-packages/fenics
RUN add-apt-repository ppa:ubuntugis/ppa

# Install libraries used by DTCCore
RUN apt-get update && apt-get install -y \
    locales \
    sudo \
	build-essential \
	cmake\
    nlohmann-json3-dev \
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
    colordiff\
    libuuid1 \
    #lftp \
    fenics \
    gdal-bin \
    libgdal-dev \ 
    git \
    #git-lfs \
    #dos2unix \
    nano \
    rsync \
    wget \
    moreutils \
    unzip \
    jq \
    uuid-dev

# Get and compile LASZip
RUN git clone -n https://github.com/LASzip/LASzip.git
RUN cd LASzip && git checkout 585a940c8d80f039fd1294ecd1411440938d7241 && ./autogen.sh && ./configure && make all -j 4 && make install
RUN mkdir /usr/local/include/laszip/ && mv /usr/local/include/las*.hpp /usr/local/include/laszip/

# Get and compile libLAS with LASZip suppport

#RUN git clone https://github.com/libLAS/libLAS.git && cd libLAS && mkdir build && cd build && cmake .. -DWITH_LASZIP=TRUE && make all -j 4 && make install
#RUN ln -s /usr/local/lib/liblas.so.3 /usr/lib && ln -s /usr/local/lib/liblaszip.so.4 /usr/lib
#Temp fix since make install has broken since 5/5/2020
RUN git clone https://github.com/libLAS/libLAS.git && cd libLAS && mkdir build && cd build && cmake .. -DWITH_LASZIP=TRUE && make all -j 4 && cp ./bin/Release/liblas.so* /usr/lib/ && cp ../include/* /usr/local/include/ -r

# Get and compile VTK 7.1

RUN wget https://www.vtk.org/files/release/7.1/VTK-7.1.1.tar.gz && tar xvzf VTK-7.1.1.tar.gz && cd VTK-7.1.1/ && mkdir build && cd build && cmake -DVTK_Group_Rendering=OFF -DVTK_BUILD_ALL_MODULES_FOR_TESTS:BOOL=OFF -DVTK_Group_StandAlone=OFF -DModule_vtkCommonCore:BOOL=ON -DModule_vtkCommonDataModel:BOOL=ON -DModule_vtkIOXML:BOOL=ON .. && make all -j 4 && make install
#-DModule_vtkFiltersGeneral:BOOL=ON

# Add /usr/local/lib into LD
RUN ldconfig /usr/local/lib

# Get and compile assimp latest (27th May 2020)

#RUN git clone https://github.com/assimp/assimp && cd assimp && git checkout 1427e67b54906419e9f83cc8625e2207fbb0fcd5 && mkdir build && cd build && cmake .. && make all -j 4 && make #install && ln -s /assimp/build/bin/libassimp.so.5 /usr/lib/x86_64-linux-gnu/


#RUN apt-get update && apt-get upgrade -y

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
#ENTRYPOINT ["/bin/bash", "-l", "-c"]
#CMD ["/bin/bash", "-i"]


FROM base as dev
USER $USER
ENTRYPOINT ["/bin/bash", "-l", "-c"]
CMD ["/bin/bash", "-i"]

FROM base as clion
USER root
########################################################
# Essential packages for remote debugging and login in
########################################################

RUN apt-get update && apt-get upgrade -y && apt-get install -y \
    apt-utils gcc g++ openssh-server cmake build-essential gdb gdbserver rsync vim

ADD . /home/dtcc/core
WORKDIR /home/dtcc/core

# Taken from - https://docs.docker.com/engine/examples/running_ssh_service/#environment-variables

#RUN mkdir /var/run/sshd
RUN echo 'root:root' | chpasswd
RUN sed -i 's/PermitRootLogin prohibit-password/PermitRootLogin yes/' /etc/ssh/sshd_config

# SSH login fix. Otherwise user is kicked off after login
RUN sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd

ENV NOTVISIBLE "in users profile"
RUN echo "export VISIBLE=now" >> /etc/profile

# 22 for ssh server. 7777 for gdb server.
EXPOSE 22 7777

RUN useradd -ms /bin/bash debugger
RUN echo 'debugger:pwd' | chpasswd

########################################################

CMD ["/usr/sbin/sshd", "-D"]


