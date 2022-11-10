# -*- coding: utf-8 -*-
from setuptools import setup, find_packages
from pathlib import Path
import os, pathlib, sys
from pybind11.setup_helpers import Pybind11Extension, build_ext


packages = [
    "pybuilder",
]

package_data = {"": ["*"]}

install_requires = [
    "numpy>=1.23.3,<2.0.0",
    "pybind11>=2.10.0,<3.0.0",
    "Fiona>=1.8.0<2.0.0",
    "Shapely>=1.8.0<2.0.0",
]

setup_kwargs = {
    "name": "pybuilder",
    "version": "0.2.8",
    "description": "python bindings for dtcc-builder",
    "author": "Dag Wästberg",
    "author_email": "dwastberg@gmail.com",
    "url": "https://gitlab.com/dtcc-platform/dtcc-builder",
    # 'packages': packages,
    "package_data": package_data,
    "install_requires": install_requires,
    "python_requires": ">=3.8,<4.0",
}


def build(setup_kwargs):
    project_root = str((Path(__file__).parent / "..").resolve())
    ext_modules = [
        Pybind11Extension(
            "_pybuilder",
            ["pybindings/src/pybind_builder.cpp"],
            include_dirs=[
                os.path.join(project_root, "include"),
                os.path.join(project_root, "external"),
                "/usr/include/nlohmann/",
                "/usr/local/include",
                "/usr/include/eigen3",
                "/usr/include",
            ],
            extra_compile_args=["-std=c++1y"],
            extra_link_args=[
                "-luuid",
                "-lshp",
                "-llas",
                "-lstdc++fs",
                "-ldolfin",
                "-ltriangle",
            ],
        ),
    ]
    setup_kwargs.update({"ext_modules": ext_modules, "zip_safe": False})


# build(setup_kwargs)

setup(**setup_kwargs, packages=find_packages(".", exclude=["pybuilder.test"]))
