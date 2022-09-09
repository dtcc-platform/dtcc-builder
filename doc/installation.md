# Installation

## Getting started

These instructions will get you up and running on your local machine
for development and testing purposes.

### Downloading the software

To download the software, clone the repository by the following command:

    git clone https://gitlab.com/dtcc-platform/dtcc-builder.git -b develop

Alternatively, you may want to use the SSH protocol:

    git clone git@gitlab.com:dtcc-platform/dtcc-builder.git

This will create a directory named `dtcc-builder` containing the full
source code.

**Note:** If you are using Windows, you might first want to make sure
that Git does not convert Unix-style file endings on checkout. This
can be accomplished by:

    git config --global dtcc-builder.autocrlf false

### Downloading demo data

To download demo data, there is an open demo dataset (public) and one
available only for development/research purposes. Enter the `data`
directory and issue the following command:

    ./dtcc-download-demo-data-public

For the typical new user, this should be done inside the
container. Obtaining the public dataset requires no credentials.

The public demo data currently contains the following two datasets:

* `Eman2021`: A residential area in the vicinity of the river Emån in
southern Sweden (Småland)
* `Majorna2021`: A residential area within the City of Gothenburg in
  Sweden

The public demo data is kindly provided by Lantmäteriet and published at

    ftp://download-open.lantmateriet.se/GEO-demoomrade

For accesing the private demo datasets you will need a set of
credentials supplied by the development team. Contact us for more
info.

### Building and installing

To build, use a standard out-of-source CMake build by issuing the
following commands from the top level directory:

    mkdir build
    cd build
    cmake ..
    make -j
	make install

This will build and install all programs into the top level `bin`
directory.

**Note:** These commands should be issued *inside* the container.
