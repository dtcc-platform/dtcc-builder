# -*- coding: utf-8 -*-
from setuptools import setup, find_packages, find_namespace_packages
from pathlib import Path
import os, pathlib, sys
from pybind11.setup_helpers import Pybind11Extension, build_ext


packages = [
    "dtccpybuilder",
]

package_data = {"": ["*.so", "*.dll", "*.dylib"]}

install_requires = [
    "numpy>=1.23.3,<2.0.0",
    "pybind11>=2.10.0,<3.0.0"
]

setup_kwargs = {
    "name": "dtcc-pybuilder",
    "version": "0.5.1",
    "description": "python bindings for dtcc-builder",
    "author": "Dag Wästberg",
    "author_email": "dwastberg@gmail.com",
    "url": "https://gitlab.com/dtcc-platform/dtcc-builder",
    # 'packages': packages,
    "packages":find_namespace_packages(include=['dtcc.*']),
    "package_data": package_data,
    "install_requires": install_requires,
    "python_requires": ">=3.8,<4.0",
}


# def build(setup_kwargs):
#     project_root = str((Path(__file__).parent / "..").resolve())
#     ext_modules = [
#         Pybind11Extension(
#             "_pybuilder",
#             ["dtcc-pybindings/src/pybind_builder.cpp"],
#             include_dirs=[
#                 os.path.join(project_root, "include"),
#                 os.path.join(project_root, "external"),
#                 "/usr/include/nlohmann/",
#                 "/usr/local/include",
#                 "/usr/include/eigen3",
#                 "/usr/include",
#             ],
#             extra_compile_args=["-std=c++1y"],
#             extra_link_args=[
#                 "-luuid",
#                 "-lshp",
#                 "-llas",
#                 "-lstdc++fs",
#                 "-ldolfin",
#                 "-ltriangle",
#             ],
#         ),
#     ]
#     setup_kwargs.update({"ext_modules": ext_modules, "zip_safe": False})


# build(setup_kwargs)

setup(**setup_kwargs)
