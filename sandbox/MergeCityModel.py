# Testing algorithm for merging polygons on Hammarkullen
# Anders Logg 2020-05-04

import json
from MergePolygons import *

minimalBuildingDistance = 0.5

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
#polygons = polygons[:100]
#polygons = polygons[4:5] + polygons[74:75]

#test0 = [24,70,228,160,313,296,167,125,230,333,52,119,187,20,101,308,154,170,85,207,195,118,326,251,257,298,210,253,77,291,297,171,245,172,81,64,54,314,106,148,129,19,290,16,275,276,196,227,272,153,206]

#polygons = [polygons[i] for i in test0]
#polygons = [polygons[523]]

# Plot original model
figure()
PlotPolygons(polygons, style='-', labels=False)
title('Original model')

# Replace self-intersecting polygons with convex hull
#for i in range(len(polygons)):
#    if not CheckPolygon(polygons[i], 0.5*minimalBuildingDistance, 0.5):
#        print('Bad polygon:', i)
#        polygons[i] = ConvexHull(polygons[i])
#        #PlotPolygons([polygons[i]], style='--', labels=False)
#        #PlotLabel(polygons[i], str(i))

# Plot cleaned model
figure()
PlotPolygons(polygons, style='-', labels=False)
title('Cleaned model')

# Create queue of indices to check
indices = list(i for i in range(len(polygons)))

#indices = []

# Used for debugging
stop = False

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
        if len(polygons[j]) == 0:
           continue;

        # Compute squared distance between polygons
        Pi = polygons[i]
        Pj = polygons[j]
        if len(Pi) == 0 or len(Pj) == 0: continue
        d = sqrt(SquaredDistancePolygonPolygon(Pi, Pj))

        #  Check if distance is smaller than the tolerance
        if d < minimalBuildingDistance:

            print('CityModelGenerator: Buildings %d and %d are too close, merging' % (i, j))

            #print(Pi)
            #print(Pj)

            # Compute merged polygon
            mergedPolygon = MergePolygons([Pi, Pj], minimalBuildingDistance)

            #polygons = [mergedPolygon]
            #PlotPolygons(polygons)
            #show()

            # Replace Pi, erase Pj and add Pi to queue
            polygons[i] = mergedPolygon
            polygons[j] = []
            indices.append(i)

        if stop: break
    if stop: break


# Extract non-empty polygons
mergedPolygons = []
for i, polygon in enumerate(polygons):
    if len(polygon) > 0:
        mergedPolygons.append(polygon)
        #print(i, ' --> ', len(mergedPolygons) - 1)

# Plot simplified model
figure()
PlotPolygons(mergedPolygons, style='-', labels=False)
title('Simplified model')

show()
