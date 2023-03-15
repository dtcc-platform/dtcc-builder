# Python API

## overview

dtcc-builder provides a python API that wraps the core C++ library that is automatically installed together with the C++ tools by cmake.  

## Example
```python
from dtccpybuilder.Parameters import load_parameters
from dtccpybuilder.CityModel import CityModel
from dtccpybuilder.PointCloud import PointCloud
from dtccpybuilder.ElevationModel import ElevationModel

# load parameters from a JSON 
parameters = load_parameters("path/to/Parameters.json") 
# load default parameter values
parameters = load_parameters()

cm = CityModel("path/to/footprints.shp",parameters)
pc = PointCloud("path/to/pointcloud.las")

pc.
```

