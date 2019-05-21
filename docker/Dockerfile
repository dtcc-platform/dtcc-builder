# This file builds the Docker image for VCCore

# Use Phusion base image (minimal Docker-friendly Ubuntu)
FROM phusion/baseimage:0.11

# Set some variables
ENV USER vcuser
ENV DIR vccore
ENV HOME /home/$USER

# Update and upgrade
RUN apt-get update && apt-get upgrade -y

# Install basic dependencies
RUN apt-get -y install locales sudo

# Install libraries used by VCCore
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
