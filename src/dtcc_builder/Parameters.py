import json
from pathlib import Path

_parameters = {}

_parameters["ModelName"] = "DTCC"
_parameters["AutoDomain"] = True
_parameters["GenerateSurfaceMeshes"] = True
_parameters["GenerateVolumeMeshes"] = True
_parameters["WriteJSON"] = True
_parameters["WriteVTK"] = True
_parameters["WriteSTL"] = True
_parameters["WriteOBJ"] = False
_parameters["WriteProtobuf"] = True

_parameters["Debug"] = False

_parameters["GroundSmoothing"] = 5

_parameters["DomainMargin"] = 10.0
_parameters["X0"] = 0.0
_parameters["Y0"] = 0.0
_parameters["XMin"] = 0.0
_parameters["YMin"] = 0.0
_parameters["XMax"] = 0.0
_parameters["YMax"] = 0.0
_parameters["ElevationModelResolution"] = 1.0
_parameters["MinBuildingDistance"] = 1.0
_parameters["MinVertexDistance"] = 1.0
_parameters["GroundMargin"] = 1.0
_parameters["MeshResolution"] = 10.0
_parameters["DomainHeight"] = 100.0
_parameters["GroundPercentile"] = 0.1
_parameters["RoofPercentile"] = 0.9
_parameters["OutlierMargin"] = 2.0
_parameters["MinBuildingSize"] = 15.0
_parameters["UUIDField"] = "uuid"
_parameters["HeightField"] = ""

_parameters["StatisticalOutlierRemover"] = True
_parameters["OutlierNeighbors"] = 5
_parameters["OutlierSTD"] = 1.5

_parameters["RANSACOutlierRemover"] = True
_parameters["RANSACOutlierMargin"] = 3.0
_parameters["RANSACIterations"] = 250

_parameters["NaiveVegitationFilter"] = True
_parameters["DataDirectory"] = ""
_parameters["BuildingsFileName"] = "PropertyMap.shp"
_parameters["PointCloudDirectory"] = ""
_parameters["OutputDirectory"] = ""


def load_parameters(file_path=None, project_path="."):
    global _parameters
    #print(file_path)
    if file_path is None:
        p = _parameters
    else:
        file_path = Path(file_path)
        if file_path.is_dir():
            file_path = file_path / "Parameters.json"
        if not file_path.exists():
            print(f"Parameters file {file_path} does not exist, using default parameters")
            p = _parameters
        else:
            with open(file_path) as src:
                loaded_parameters = json.load(src)
            _parameters.update(loaded_parameters)
            p = set_directories(_parameters, project_path)
    return p

def set_directories(p, project_path):
    if p["DataDirectory"] == "":
        p["DataDirectory"] = project_path
    if p["OutputDirectory"] == "":
        p["OutputDirectory"] = p["DataDirectory"]

    p["DataDirectory"] = Path(p["DataDirectory"])
    p["OutputDirectory"] = Path(p["OutputDirectory"])
    p["OutputDirectory"].mkdir(parents=True, exist_ok=True)
    if p["PointCloudDirectory"] == "":
        p["PointCloudDirectory"] = p["DataDirectory"]
    else:
        p["PointCloudDirectory"] = Path(p["PointCloudDirectory"])
        if not p["PointCloudDirectory"].is_absolute():
            p["PointCloudDirectory"] = p["DataDirectory"] / \
                p["PointCloudDirectory"]

    return p