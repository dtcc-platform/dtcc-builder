# Edit the city model for the simulation demo before generating
# the mesh to make the city model more interesting by extending
# the heights of some of the buildings.

import json

new_heights = {8: 40.0, 21: 60.0}

print("Adjusting building heights for simulation demo...")

path = "../data/HelsingborgResidential2022/CityModel.json"
with open(path, "r") as f:
    data = json.load(f)

for i, h in new_heights.items():
    b = data["Buildings"][i]
    data["Buildings"][i]["Height"] = h
    for j in range(len(b["RoofPoints"])):
        data["Buildings"][i]["RoofPoints"][j]["z"] = b["GroundHeight"] + h
    # print(b["GroundHeight"])
    # print(b["Height"])
    # print([p["z"] for p in b["RoofPoints"]])
    # print("")

with open(path, "w") as f:
    json.dump(data, f, indent=4)
