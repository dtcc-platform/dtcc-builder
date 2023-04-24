import dtcc_io as io
import dtcc_builder.build

cm = io.load_citymodel(
    "/home/cit-wastberg/dtcc-platform/dtcc-builder/data/HelsingborgResidential2022/PropertyMap.shp"
)
pc = io.load_pointcloud(
    "/home/cit-wastberg/dtcc-platform/dtcc-builder/data/HelsingborgResidential2022/PointCloud.las"
)

cm = dtcc_builder.build.extract_buildingpoints(cm, pc)
