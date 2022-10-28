import _pybuilder
import os

def ReadLasFiles(las_path, extra_data = True):
    las_path = str(las_path)
    
    print(f"loading las from {las_path}")
    if os.path.isdir(las_path):
        print("\tloading as dir")
        pc = _pybuilder.LASReadDirectory(las_path, extra_data)
    elif os.path.isfile(las_path):
        print("\tloading as file")
        pc = _pybuilder.LASReadFile(las_path,extra_data)
    return pc

def GetLasBounds(las_path):
    pass
 