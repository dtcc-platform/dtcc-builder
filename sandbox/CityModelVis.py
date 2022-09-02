import json
from MergePolygons import *

minimalBuildingDistance = 0.5
PLOT_LABELS = False
PLOT_ALGORITHM = True

def extractPolygons(data):
    polygons=[]
    for building in data['Buildings']:
        polygon = []
        for p in building['Footprint']:
            polygon.append((p['x'], p['y']))
        polygon = array(polygon)
        polygons.append(polygon)
    return polygons

# Read building footprints for Hammarkullen
#with open('../data/Hammarkullen/CityModel.json') as f:
#with open('..data/Lund2022/CityModel.json') as f:
with open('E:\GitHub_Projects\dtcc-builder\data\Lund2022\CityModel.json') as f:
    data = json.loads(f.read())
polygons = extractPolygons(data)

with open('E:\GitHub_Projects\dtcc-builder\data\Lund2022\CityModelSimple.json') as f_simple:
    data_simple=json.loads(f_simple.read())
polygons_simple = extractPolygons(data_simple)

figure()
PlotPolygons(polygons, style='-', labels=PLOT_LABELS)
PlotPolygons(polygons_simple,style='-', labels=PLOT_LABELS)
title('Overlapping models')

figure()
PlotPolygons(polygons, style='-', labels=PLOT_LABELS)
title('Original model')
#show()

figure()
PlotPolygons(polygons_simple,style='-', labels=PLOT_LABELS)
title('Simple model')
show()

# Plot original model
#figure()
#PlotPolygons(polygons, style='-', labels=PLOT_LABELS)
#title('Original model')

#with open('E:\GitHub_Projects\dtcc-builder\data\Lund2022\CityModelSimple.json') as f_simple:
#    data_simple=json.loads(f_simple.read())
#
#polygons_simple = extractPolygons(data_simple)
#figure()
#PlotPolygons(polygons_simple,style='-', labels=PLOT_LABELS)
#show()
