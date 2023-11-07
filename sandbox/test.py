import dtcc_model as model
import dtcc_io as io
import dtcc_builder as builder



def main():
    # Get parameters
    p = builder.parameters.default()

    p["auto_domain"] = True
    #[102000.0, 102999.99] x [6213004.15744457, 6213999.99]
    p["x_min"] = 102000.0
    p["x_max"] = 102000.0 + 100

    p["y_min"] = 6213004.15744457
    p["y_max"] = 6213004.15744457 + 100
    buildings_path = "/home/auth-georspai/dtcc-builder/data/HelsingborgResidential2022/PropertyMap.shp"
    pointcloud_path = "/home/auth-georspai/dtcc-builder/data/HelsingborgResidential2022"

    origin, bounds = builder.calculate_bounds(buildings_path, pointcloud_path, p)
    print(bounds)

       # Load city
    city = io.load_city(
        buildings_path,
        uuid_field=p["uuid_field"],
        height_field=p["height_field"],
        bounds=bounds,
    )

    # Load point cloud
    point_cloud = io.load_pointcloud(pointcloud_path, bounds=bounds)

    # Build city
    city = builder.build_city(city, point_cloud, bounds, p)




if __name__== "__main__":
    main()