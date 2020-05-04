# Testing algorithm for merging polygons on Hammarkullen
# Anders Logg 2020-05-04

import json
from MergePolygons import *

minimalBuildingDistance = 10.0

# Read building footprints for Hammarkullen
with open('../data/Hammarkullen/CityModel.json') as f:
    data = json.loads(f.read())

# Extract polygons
polygons = []
for building in data['Buildings']:
    polygon = []
    for p in building['Footprint']:
        polygon.append((p['x'], p['y']))
    polygon = array(polygon)
    polygons.append(polygon)

# FIXME: Testing
polygons = polygons[:10]


# Create queue of indices to check
indices = list(i for i in range(len(polygons)))

# Process queue until empty
while len(indices) > 0:

    # Pop index of next polygon to check
    i = indices.pop(0)

    # Iterate over all other polygons
    for j in range(len(polygons)):

        # Skip polygon itself
        if i == j:
            continue

        # Skip if polygon has zero size (merged with other polygon)
        if len(polygons) == 0:
           continue;

        # Compute squared distance between polygons
        Pi = polygons[i]
        Pj = polygons[j]
        d = sqrt(SquaredDistancePolygonPolygon(Pi, Pj))

        #  Check if distance is smaller than the tolerance
        if d < minimalBuildingDistance:

            print('CityModelGenerator: Buildings %d and %d are too close, merging' % (i, j))

            print(Pi)
            print(Pj)

            # Compute merged polygon
            mergedPolygon = MergePolygons([Pi, Pj], minimalBuildingDistance)

            # Replace Pi, erase Pj and add Pi to queue
            polygons[i] = mergedPolygon
            polygons[j] = []
            indices.append(i)
