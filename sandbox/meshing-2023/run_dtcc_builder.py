#!/usr/bin/env python3


import json, subprocess, glob, shutil, os, datetime, copy

binary_dtcc_builder = "dtcc-build"

# Parameters for mesh resolution study
parameter_name = 'mesh_resolution'
parameter_value = [ 1.0 , 2.0, 4.0 ]
# parameter_value = [ 10.0, 16.0 , 32.0 ]
parameter_value.sort(reverse=True)
parameter_file = "./parameters.json"
if not os.path.isfile(parameter_file):
    print(f'Invalid parameters file: {parameter_file}')
    exit()

output_directory = "./output"


if not os.path.isdir(output_directory):
    print(f'Invalid output directory: {output_directory}')
    exit()

# Read parameters
with open(parameter_file,'r') as f:
    parameters = json.load(f)


for resolution in parameter_value:
    print(f'\033[1;32m Running dtcc_builder with {parameter_name}: {resolution}\033[0;37m')

    #Update parameters
    parameters[parameter_name] = resolution
    parameters["output_directory"] = os.path.join(output_directory, "mesh_" + str(resolution).replace(".","_"))
    with open(parameter_file,'w') as f:
        json.dump(parameters, f , indent= 2)

    if not os.path.exists(parameters['output_directory']):
        os.mkdir(parameters['output_directory'])

    process = subprocess.Popen([binary_dtcc_builder, parameter_file])
    process.wait()


