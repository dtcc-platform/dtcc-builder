import time

import dtcc_io as io
from dtcc_model.gridfield import GridField2D
from dtcc_builder import _pybuilder, builder_datamodel
import numpy as np
import affine
import rasterio

pc = io.load_pointcloud(
    "/home/cit-wastberg/dtcc-platform/dtcc-builder/data/HelsingborgResidential2022/PointCloud.las"
)

cell_size = 0.5
start_time = time.time()
b_pc = builder_datamodel.create_builder_pointcloud(pc)
b_gf = _pybuilder.GenerateElevationModel(b_pc, cell_size, [2, 9])
b_gf = _pybuilder.SmoothElevation(b_gf, 5)
print("Time to generate elevation model: ", time.time() - start_time)

gf = GridField2D()
shape = b_gf.Grid.xsize, b_gf.Grid.ysize
print(shape)
gf.grid = np.array(b_gf.values).reshape(shape)
gf.transform = rasterio.transform.from_origin(
    pc.bounds[3], pc.bounds[0], cell_size, cell_size
)
io.save_gridfield(gf, "test5.tif")
