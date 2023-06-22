import dtcc_io as io
import dtcc_builder.build

city = io.load_city(
    "/home/cit-wastberg/dtcc-platform/dtcc-builder/data/HelsingborgResidential2022/PropertyMap.shp"
)
pc = io.load_pointcloud(
    "/home/cit-wastberg/dtcc-platform/dtcc-builder/data/HelsingborgResidential2022/PointCloud.las"
)

city = dtcc_builder.build.extract_buildingpoints(city, pc)
