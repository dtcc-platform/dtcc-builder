# Installation

## Downloading the software

To download the software, clone the Git repository by the following command:

    git clone https://gitlab.com/dtcc-platform/dtcc-builder.git

Alternatively, you may want to use the SSH protocol:

    git clone git@gitlab.com:dtcc-platform/dtcc-builder.git

This will create a directory named `dtcc-builder` containing the full
source code of DTCC Builder.

## Installing dependencies

DTCC Builder depends on a number of open-source libraries. The easiest
way to install these dependencies is to use the provided Docker image
for DTCC Builder. To build and start the DTCC Docker image
(container), enter the `docker` directory and issue the following two
commands:

    bin/dtcc-build
    bin/dtcc-start

The first of these two commands will build a Docker image and
container for DTCC Builder and the second command will start and
attach to the container.

## Building the software

DTCC Builder is implemented in C++ and uses
[CMake](https://cmake.org/) for configuration and build. To create a
standard out-of-source CMake build, issue the following commands from
the top-level directory:

    mkdir build
    cd build
    cmake ..
    make -j
	make install

This will build and install all programs into the top-level `bin`
directory.

> **Note:** These commands should be run from inside the Docker
> container provided for DTCC Builder.

## Downloading demo data

DTCC builder requires data to run. The demo data are not part of the
repository and must be downloaded separately. To download public demo
data, enter the `data` directory and issue the following command:

    ./dtcc-download-demo-data-public

> **Note:** The command `dtcc-download-demo-data` provides additional
> datasets that may not be shared publicly (because of license
> restrictions from the data owners).
