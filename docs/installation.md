# Installation

DTCC has an external dependency on the libgeos library. This library must be install
using your package manager on Linux, or homebrew on MacOS.

DTCC Builder can be easily installed using [`pip`](https://pypi.org/project/pip/).


To install from the Python Package Index (PyPI):

    pip install dtcc-builder

To install from the source directory:

    pip install .

> **NOTE**: Fix for conflicting versions of `pip` and `python`
>
> Sometimes `pip` and `python` may be out of sync which means that `pip` will
> install a package in a location where it will not be found by `python`.
> It is therefore safer to replace the `pip` command by `python -m pip`:
>
>     python -m pip install [ package-name or . ]

> **NOTE**: Fix for broken setuptools in Ubuntu 22.04
>
> A bug in Ubuntu 22.04 prevents [PEP621](https://peps.python.org/pep-0621/)
> compliant Python projects from installing properly with `pip`, resulting in
> package name and version number `UNKNOWN-0.0.0`.
> To fix this, run the following commmand before `pip install`:
>
>     export DEB_PYTHON_INSTALL_LAYOUT=deb_system

## Docker image

DTCC Builder depends on a number of open-source libraries. The easiest
way to install these dependencies is to use the provided Docker image
for DTCC Builder. To build and start the DTCC Docker image
(container), enter the `docker` directory and issue the following two
commands:

    ./docker-build-container
    ./docker-start-container

The first of these two commands will build a Docker image and
container for DTCC Builder and the second command will start
the container.

## Demo data

DTCC builder requires data to run. The demo data are not part of the
repository and must be downloaded separately. To download public demo
data, enter the `data` directory and issue the following command:

    ./download-demo-data

> **Note:** The command `download-demo-data=private` provides additional
> datasets that may not be shared publicly (because of license
> restrictions from the data owners).
