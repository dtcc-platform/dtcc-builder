import json
import numpy as np
from pathlib import Path
from shapely.geometry import LineString, Polygon, MultiPoint, Point
import math

from roof_segmentation import get_roof_segments

def get_roof_data(building):
    footprint = [ [f['x'], f['y']] for f in building['Footprint']]
    footprint = Polygon(footprint)
    points = [[p['x'], p['y'],p['z']] for p in building['RoofPoints']]
    points = np.array(points)
    return (footprint, points)

def segment_cityModel_roofs(cityModel):
    buildings = cityModel['Buildings']
    for building in buildings:
        footprint, roofpoints = get_roof_data(building)
        if len(roofpoints) > 20 :
            roof_segments = get_roof_segments(roofpoints,8,)
            if len(roof_segments) > 0:
                pass


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        #print("Usage: dtcc-segment-roof.py <parameters.json>")
        #sys.exit(1)
        parameters = "demo/Majorna2021.json"
    else:
        parameters = sys.argv[1]
    with open(parameters) as src:
        params = json.load(src)
        data_dir = Path(params["DataDirectory"])
    if not data_dir.is_dir():
        print("Data directory does not exist: {}".format(data_dir))
        sys.exit(1)
    cityModel = data_dir / "CityModel.json"
    if not cityModel.is_file():
        print("City model does not exist: {}".format(cityModel))
        sys.exit(1)
    with open(cityModel) as src:
        city = json.load(src)
    segmented_city_model = segment_cityModel_roofs(city)
    


