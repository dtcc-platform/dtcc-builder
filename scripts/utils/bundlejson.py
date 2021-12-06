#!/bin/python3
# Simple little script to bundle JSON data for the online
# viewer for Andreas Rudena for the Centre Day demo.
# Anders Logg 2021-12-05

import json

# Load city model
print('Loading city model...')
with open('CityModel.json') as f:
    cityModel = json.load(f)

# Load ground surface
print('Loading ground surface...')
with open('GroundSurface.json') as f:
    groundSurface = json.load(f)

# Remove building points
print('Removing building points...')
for building in cityModel['Buildings']:
    del building['GroundPoints']
    del building['RoofPoints']

# Add ground surface
cityModel['GroundSurface'] = groundSurface

# Write to file
print('Writing updated city model...')
with open('CityModelSpecial.json', 'w') as f:
    json.dump(cityModel, f)
