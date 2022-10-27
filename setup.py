# -*- coding: utf-8 -*-
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext


packages = \
['pybuilder']

package_data = \
{'': ['*'], 'pybuilder': ['include/*', 'pybind/src/*']}

install_requires = \
['numpy>=1.23.3,<2.0.0', 'pybind11>=2.10.0,<3.0.0']

setup_kwargs = {
    'name': 'pybuilder',
    'version': '0.1.0',
    'description': 'python bindings for dtcc-builder',
    'author': 'Dag Wästberg',
    'author_email': 'dwastberg@gmail.com',
    'url': 'https://gitlab.com/dtcc-platform/dtcc-builder',
    'packages': packages,
    'package_data': package_data,
    'install_requires': install_requires,
    'python_requires': '>=3.8,<4.0',
}


def build(setup_kwargs):
    ext_modules = [
        Pybind11Extension("_pybuildier", ["pybind/src/pybind_builder.cpp"],
        include_dirs = ["include","pybind/include"],
        extra_compile_args=['-std=c++1y']),
    ]
    setup_kwargs.update({
        "ext_modules": ext_modules,
        "zip_safe": False
    })

build(setup_kwargs)

setup(**setup_kwargs)