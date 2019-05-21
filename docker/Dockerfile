# This file builds the Docker image for VCCore

# Use Phusion base image (minimal Docker-friendly Ubuntu)
FROM phusion/baseimage:0.11

# Set some variables
ENV USER vcuser
ENV HOME /home/$USER

# Install basic dependencies
#RUN apt-get update && apt-get install -y \
#    nlohmann-json-dev \
#    libshp-dev \
#    liblas-dev \
#    liblas-c-dev \
#    libpugixml-dev \
#    libproj-dev \
#    libtriangle-dev

# Install FEniCS
#RUN add-apt-repository ppa:fenics-packages/fenics
#RUN apt-get update && apt-get install -y fenics

# Add user and change to user
RUN useradd -m -s /bin/bash -G sudo $USER
USER $USER
WORKDIR $HOME

# Generate welcome message printed at login
COPY Welcome $HOME/.welcome
RUN echo "cat $HOME/.welcome" >> $HOME/.bashrc

ENTRYPOINT ["bash", "-l"]
