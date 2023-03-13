# This is the base Dockerfile for DTCC Platform. This file is part
# of the dtcc-model repository (which resides as a submodule in all
# other repositories). Docker images for all other repositories are
# built on top of this base image by adding specific dependencies
# via the install_dependencies.sh script.

# Use Phusion base image (minimal Docker-friendly Ubuntu)
FROM phusion/baseimage:jammy-1.0.1 as base

# Set some variables
ENV USER dtcc
ENV HOME /home/$USER

# Install libraries used by other repos
COPY install_scripts/ /
#COPY docker/createbanner.sh /

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y locales  sudo dialog apt-utils
RUN TZ=Etc/UTC apt-get -y install tzdata

RUN ./install_dependencies.sh

# Add /usr/local/lib into LD
RUN ldconfig /usr/local/lib

# Add user and change to user
RUN useradd -m $USER -G sudo && \
    echo "$USER:$USER" | chpasswd && \
    echo "$USER ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers
ARG _uid
ARG gid
RUN usermod -u ${_uid} dtcc
RUN groupmod -g ${gid} dtcc

#RUN touch banner
#RUN ./createbanner.sh
#USER $USER

# Generate welcome message printed at login
#RUN touch banner
#RUN ./createbanner.sh
#COPY banner /
#COPY banner $HOME/.welcome
#COPY *docker/Welcome $HOME/.welcome
#RUN echo "cat $HOME/.welcome" >> $HOME/.bashrc

# Core FastAPI port
EXPOSE 8090

RUN mkdir /shared_data
RUN chown $USER:$USER /shared_data
RUN chown $USER:$USER /home/$USER

# Image setup
FROM base as dev
USER $USER
ENTRYPOINT ["/bin/bash", "-l", "-c"]
CMD ["/bin/bash", "-i"]
