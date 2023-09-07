# Copyright (C) 2022 Dag Wästberg
# Licensed under the MIT License
#
# Modified by Anders Logg 2023
#
# This is the main command line script provided by DTCC Builder.
# It builds city models and/or meshes from raw data in a a given
# project directory.
#
# This script is responsible for parsing command line arguments
# and handle logic around data directories and parameters. The
# actual building is delegated to the main build() function.

import sys, os, re, json, argparse
from pathlib import Path


import dtcc_io as io
import dtcc_builder as builder
from dtcc_builder.logging import info, warning

PARAMETERS_FILE = "parameters.json"


def parse_command_line():
    """
    Parse command line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command line arguments.
    """

    # Create parser
    parser = argparse.ArgumentParser(
        prog="dtcc-builder",
        description="Build LOD1 city mesh(es) from building footprints and pointcloud",
    )

    # Add main arguments
    parser.add_argument("path", nargs="?", default=os.getcwd())
    # parser.add_argument("--city-only", action="store_true")
    # parser.add_argument("--mesh-only", action="store_true")

    # Add parameter arguments (note special handling of booleans)
    for key, value in builder.parameters.default().items():
        if isinstance(value, bool):
            parser.add_argument(
                f"--{key}",
                dest=key,
                action="store_true",
                default=None,
                help=f"Turn on {key}",
            )
            parser.add_argument(
                f"--no-{key}",
                dest=key,
                default=None,
                action="store_false",
                help=f"Turn off {key}",
            )
        else:
            parser.add_argument(f"--{key}", default=None, type=type(value))

    return parser.parse_args()


def load_parameters(args):
    """
    Load parameters and override by command line arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command line arguments.

    Returns
    -------
    dict
        Loaded and overridden parameters.
    """

    # Load parameters from or directory if provided, else use default parameters
    parameters = builder.parameters.default()
    if args.path is None:
        info("Using default parameters")
    else:
        path = Path(args.path)
        if path.is_dir():
            path = path / PARAMETERS_FILE
        if path.exists():
            info(f"Loading parameters from {path}")
            with open(path) as f:
                _parameters = json.load(f)
                parameters.update(_parameters)
        else:
            warning(f"Unable to load {path}; using default parameters")

    # Override parameters with command line arguments
    for key, value in parameters.items():
        parser_value = getattr(args, key)
        if parser_value is not None:
            info(f"Overriding parameter {key} with argument {parser_value}")
            parameters[key] = parser_value

    # Set parameters for directories
    set_directory_parameters(parameters, args.path)

    # Pretty-print parameters
    info("Printing parameters")
    keys = sorted(parameters.keys())
    n = max([len(key) for key in keys]) if keys else 0
    for key in sorted(keys):
        print(f"  {key}: {' '*(n - len(key) - 1)} {parameters[key]}")

    return parameters


def set_directory_parameters(parameters, path):
    """
    Set parameters for directories.

    Parameters
    ----------
    parameters : dict
        Dictionary of parameters to be set.
    path : str
        Path to the data directory.

    Returns
    -------
    dict
        Updated parameters dictionary.
    """
    # Shortcut
    p = parameters

    # Set data_directory
    data_directory = p["data_directory"]
    if data_directory == "":
        data_directory = path
    data_directory = Path(data_directory)

    # Set pointcloud_directory
    pointcloud_directory = p["pointcloud_directory"]
    if pointcloud_directory == "":
        pointcloud_directory = data_directory
    pointcloud_directory = Path(pointcloud_directory)

    # Set output_directory
    output_directory = p["output_directory"]
    if output_directory == "":
        output_directory = data_directory
    output_directory = Path(output_directory)

    # Set parameters
    p["data_directory"] = data_directory
    p["pointcloud_directory"] = pointcloud_directory
    p["output_directory"] = output_directory

    return p


def main():
    # Parse command line
    args = parse_command_line()

    # Load parameters
    parameters = load_parameters(args)

    # Call main build function
    builder.build(parameters)
