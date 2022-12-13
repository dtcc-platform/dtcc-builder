
import math
from roadnetwork_pb2 import RoadNetwork, Road, Vector2D, Georeference;

def GetExistingVertexIndexZoned(gVertices, newVertex, zoneIndex, zoneVertexDict):
    tol = 0.001
    if zoneIndex in zoneVertexDict:
        for vIndex in zoneVertexDict[zoneIndex]:
            v = gVertices[vIndex]
            d = Distance(v, newVertex)
            if d < tol:
                return vIndex
    return -1    

def Distance(v1, v2):
    return math.sqrt(math.pow((v1.x - v2.x),2) + math.pow((v1.y - v2.y),2)) 

def ShouldExistingIndexBeAdded(rVertices, testIndex):
    if len(rVertices) > 0:
        lastIndex = rVertices[-1]
        if(lastIndex != testIndex):
            return True
        else:
            return False    
    else:
        return True     

def GetRoadType(strType):    
    if(strType == 'VÄGKV.M'):
        roadType = 2
    elif (strType == 'VÄGBN.M'):
        roadType = 3
    elif (strType == 'VÄGGG.M'):
        roadType = 4
    elif (strType == 'VÄGBNU.M'):
        roadType = 5
    else:
        roadType = 0
    
    return roadType

def AppendNewOrExistingVertexZoned(gVertices, gVertexIndices, gVertexIndex, rVertexIndices, v, zoneIndex, zoneVertexDict):
    existingVertexindex = GetExistingVertexIndexZoned(gVertices, v, zoneIndex, zoneVertexDict)
    
    if existingVertexindex == -1: #The new vertex does not exist already
        gVertices.append(v)
        rVertexIndices.append(gVertexIndex)
        gVertexIndices.append(gVertexIndex)
        if zoneIndex in zoneVertexDict:
            zoneVertexDict[zoneIndex].append(gVertexIndex)
        else:
            zoneVertexDict[zoneIndex] = []
            zoneVertexDict[zoneIndex].append(gVertexIndex)

        gVertexIndex = gVertexIndex + 1
    else:                         #The vertex already existed in the global list 
        if(ShouldExistingIndexBeAdded(rVertexIndices, existingVertexindex)):  #Check to avoid having two of the same index in a row       
            rVertexIndices.append(existingVertexindex)
    
    return gVertexIndex        


def GetVertexZoneIndex(xMin, yMin, nx, ny, dx, dy, v):
    zoneIndex = math.floor((v.x - xMin) / dx) * ny + math.floor((v.y - yMin) / dy )
    return zoneIndex

def PreprocessZones(shpFile, xDom, yDom, nn, dd, nx):
    bb = shpFile.bounds
    bbPadding = 50.0                                            #Adding some padding to make bb slightly larger    
    xDom.extend([bb[0]-bbPadding, bb[2]+bbPadding])             #domain start at zero
    yDom.extend([bb[1]-bbPadding, bb[3]+bbPadding])             #for both directions 
    dx = (xDom[1] - xDom[0]) / nx 
    ny = math.ceil((yDom[1] - yDom[0]) / dx)
    dy = (yDom[1] - yDom[0]) / ny
    dd.extend([dx, dy])
    nn.extend([nx, ny])


def PrintResults(shouldPrint, roadNetwork, zoneVertexDict):

    if shouldPrint:
        vertexCount = 0
        with open("export/vertexExport.txt", "w") as f:
            for v in roadNetwork.vertex:
                f.write(str(v))

        with open("export/roadExport.txt", "w") as f:
            for r in roadNetwork.roads:
                f.write(str(r.vertices)+'\n')

        with open("export/zoneDict.txt", "w") as f:
            for key in zoneVertexDict:
                f.write('Zone index:' + str(key) + ' ' + str(zoneVertexDict[key])+'\n')
                vertexCount = vertexCount + len(zoneVertexDict[key])   
        print(vertexCount)    




