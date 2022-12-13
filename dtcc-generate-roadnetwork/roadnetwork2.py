import fiona
import math
import enum
import os
import time
import roadNetworkUtils as utils
from roadnetwork_pb2 import RoadNetwork, Road, Vector2D, Georeference;

os.system('clear')

print("-------- Road network reader started -------")
time1 = time.perf_counter()

#fileName = 'data/TV_subsection.shp'      #Trafikverket data
fileName = 'data/TrafikverketGbgC.shp'      #Trafikverket data
#fileName = 'data/vl_riks.shp'           #Lantmäteriet data

shpFile = fiona.open(fileName, 'r')

featureCount = sum(1 for _ in shpFile)
print("Number of features in shape file: " + str(featureCount))
nx = 10 * math.ceil(math.sqrt(featureCount))           #A good way to chose nx?

xDom, yDom, dd, nn = [], [], [], []
utils.PreprocessZones(shpFile, xDom, yDom, nn, dd, nx)
print("Number of zones: " + "[" + str(nn[0]) + "," + str(nn[1]) + " = " + str(nn[0] * nn[1]) + "]" )
zoneVertexDict = {}

mod = 100
if featureCount < 10000 and featureCount > 1000:
    mod = 1000
elif featureCount > 10000:
    mod = 10000     

exportData = True
counter = 0
roads = []
gVertices = []
gVertexIndices = []
gVertexIndex = 0

for shpObjKey, shpObject in shpFile.items(): #Looping over roads 

    roadName = "Undefined"
    roadType = 0
    geomType = shpObject['geometry']['type']
    
    #Find properties
    if(geomType == 'LineString' or geomType == 'MultiLineString'):        
        for prop in shpObject['properties'].items():
            if prop[0] == 'DETALJTYP':
                roadType = utils.GetRoadType(prop[1])     
            elif prop[0] == 'NAMN1':
                if isinstance(prop[1], str):
                    roadName = prop[1]   

    if(geomType == 'LineString'):
        rVertexIndices = []
        road = Road()
        for obj in shpObject['geometry']['coordinates']:
            v = Vector2D()
            v.x = obj[0]
            v.y = obj[1]
            zoneIndex = utils.GetVertexZoneIndex(xDom[0], yDom[0], nn[0], nn[1], dd[0], dd[1], v)
            gVertexIndex = utils.AppendNewOrExistingVertexZoned(gVertices, gVertexIndices, gVertexIndex, rVertexIndices, v, zoneIndex, zoneVertexDict)
             
            #Append to current road    
            road.vertices.extend(rVertexIndices)
            road.type = roadType
            road.name = roadName
            roads.append(road)           

    elif (geomType == 'MultiLineString'):            
        for obj in shpObject['geometry']['coordinates']:
            rVertexIndices = []
            road = Road()
            for c in obj:
                v = Vector2D()
                v.x = c[0]
                v.y = c[1]
                zoneIndex = utils.GetVertexZoneIndex(xDom[0], yDom[0], nn[0], nn[1], dd[0], dd[1], v)
                gVertexIndex = utils.AppendNewOrExistingVertexZoned(gVertices, gVertexIndices, gVertexIndex, rVertexIndices, v, zoneIndex, zoneVertexDict)
              
            #Create a road for each linestring in the multilinestring 
            road.vertices.extend(rVertexIndices)
            road.type = roadType
            road.name = roadName
            roads.append(road)

    #End lineString multiString..
#End road loop

time2 = time.perf_counter()
print("Time elapsed: "  + str(round(time2  - time1,2)))

#geoRef = Georeference()
#geoRef.crs = 'test'
#geoRef.epsg = 1
#geoRef.x0 = 0.4
#geoRef.y0 = 0.4

#Create a road network
roadNetwork = RoadNetwork()
roadNetwork.vertex.extend(gVertices)
roadNetwork.roads.extend(roads)
#roadNetwork.georeference = geoRef

#utils.PrintResults(exportData, roadNetwork, zoneVertexDict)




