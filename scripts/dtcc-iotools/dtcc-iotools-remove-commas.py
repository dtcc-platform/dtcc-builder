#!/usr/bin/env python
## VirtualCity@Chalmers: vc-iotools-removecommas
## Vasilis Naserentin 2019
## Removes commas from a csv file using fileinput. Needed for gdaltransform
import os.path
import fileinput
import sys

if len(sys.argv) < 3 or len(sys.argv) > 3:
    print("Usage: vc-removecsv filein fileout")
    exit()
filein = sys.argv[1]
fileout = sys.argv[2]
if os.path.exists(filein):
    f = open(fileout, "w")
    for line in fileinput.input(filein, inplace=False):
        string = line.replace(",", " ")
        f.write(string)
else:
    print("File " + filein + " not found")
