import json
from pathlib import Path

_parameters = {}

_parameters["model-name"] = "DTCC"
_parameters["autodomain"] = True
_parameters["generate-surface-mesh"] = True
_parameters["generate-volume-mesh"] = True
_parameters["write-json"] = True
_parameters["write-vtk"] = True
_parameters["write-stl"] = True
_parameters["write-obj"] = False
_parameters["write-protobuf"] = True

_parameters["debug"] = False

_parameters["ground-smoothing"] = 5

_parameters["domain-margin"] = 10.0
_parameters["x0"] = 0.0
_parameters["y0"] = 0.0
_parameters["x-min"] = 0.0
_parameters["y-min"] = 0.0
_parameters["x-max"] = 0.0
_parameters["y-max"] = 0.0
_parameters["elevation-model-resolution"] = 1.0
_parameters["min-building-distance"] = 1.0
_parameters["min-building-height"] = 2.5
_parameters["min-vertex-distance"] = 1.0
_parameters["ground-margin"] = 1.0
_parameters["mesh-resolution"] = 10.0
_parameters["domain-height"] = 100.0
_parameters["groun-percentile"] = 0.5
_parameters["roof-percentile"] = 0.9
_parameters["outlier-margin"] = 2.0
_parameters["min-building-size"] = 15.0
_parameters["uuid-field"] = "uuid"
_parameters["height-field"] = ""

_parameters["statistical-outlier-remover"] = True
_parameters["outlier-neighbors"] = 5
_parameters["roof-outlier-margin"] = 1.5

_parameters["ransac-outlier-remover"] = True
_parameters["ransac-outlier-margin"] = 3.0
_parameters["ransac-iterations"] = 250

_parameters["naive-vegitation-filter"] = True
_parameters["data-directory"] = ""
_parameters["buildings-filename"] = "PropertyMap.shp"
_parameters["pointcloud-directory"] = ""
_parameters["output-directory"] = ""


def load_parameters(file_path=None, project_path="."):
    global _parameters

    if file_path is None:
        p = _parameters.copy()

    else:
        file_path = Path(file_path)
        print(file_path)
        if file_path.is_dir():
            file_path = file_path / "Parameters.json"
        if not file_path.exists():
            print(
                f"Parameters file {file_path} does not exist, using default parameters"
            )
            p = _parameters.copy()
        else:
            with open(file_path) as src:
                loaded_parameters = json.load(src)
            p = _parameters.copy()
            p.update(loaded_parameters)
    p = set_directories(p, project_path)
    return p


def set_directories(p, project_path):
    if p["data-directory"] == "":
        p["data-directory"] = project_path
    if p["output-directory"] == "":
        p["output-directory"] = p["data-directory"]

    p["data-directory"] = Path(p["data-directory"])
    p["output-directory"] = Path(p["output-directory"])
    p["output-directory"].mkdir(parents=True, exist_ok=True)
    if p["pointcloud-directory"] == "":
        p["pointcloud-directory"] = p["data-directory"]
    else:
        p["pointcloud-directory"] = Path(p["pointcloud-directory"])
        if not p["pointcloud-directory"].is_absolute():
            p["pointcloud-directory"] = p["data-directory"] / p["pointcloud-directory"]

    return p
