# Installation

## Downloading the software

To download the software, clone the Git repository by the following command:

    git clone https://gitlab.com/dtcc-platform/dtcc-builder.git -b develop

Alternatively, you may want to use the SSH protocol:

    git clone git@gitlab.com:dtcc-platform/dtcc-builder.git -b develop

This will create a directory named `dtcc-builder` containing the full
source code of DTCC Builder.

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

> **Note:** DTCC Builder depends on a number of open-source
libraries. To ensure that you have all the dependencies needed, it is
recommended to work from inside the [DTCC
Docker](https://gitlab.com/dtcc-platform/dtcc-docker) container for
DTCC Builder. This documentation assumes that all commands are run
from within the Docker container.

## Downloading demo data

DTCC builder requires data to run. The demo data are not part of the
repository and must be downloaded separately. To download public demo
data, enter the `data` directory and issue the following command:

    ./dtcc-download-demo-data-public

> **Note:** The command `dtcc-download-demo-data` provides additional
datasets that may not be shared publicly (because of license
restrictions from the data owners). Access to no-public data is
provided to developers on request. Contact the maintainers of DTCC
Builder for more info.
