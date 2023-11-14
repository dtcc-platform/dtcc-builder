_parameters = {}

_parameters["model_name"] = "DTCC"
_parameters["auto_domain"] = True
_parameters["debug"] = False

_parameters["build_mesh"] = True
_parameters["build_volume_mesh"] = True

_parameters["save_protobuf"] = True
_parameters["save_json"] = False
_parameters["save_shp"] = False
_parameters["save_vtk"] = False
_parameters["save_stl"] = False
_parameters["save_obj"] = False

_parameters["ground_smoothing"] = 3

_parameters["domain_margin"] = 10.0
_parameters["x0"] = 0.0
_parameters["y0"] = 0.0
_parameters["x_min"] = 0.0
_parameters["y_min"] = 0.0
_parameters["x_max"] = 0.0
_parameters["y_max"] = 0.0
_parameters["elevation_model_resolution"] = 1.0
_parameters["elevation_model_window_size"] = 3
_parameters["min_building_detail"] = 1.0
_parameters["min_building_angle"] = 10.0
_parameters["min_building_height"] = 2.5  # sets if too small (rethink?)
_parameters["min_building_area"] = 15.0  # removes if too small (rethink?)
_parameters["min_mesh_angle"] = 30.0  # not guaranteed to converge
_parameters["max_mesh_size"] = 10.0
_parameters["ground_margin"] = 1.0
_parameters["domain_height"] = 100.0
_parameters["ground_percentile"] = 0.5
_parameters["roof_percentile"] = 0.9
_parameters["outlier_margin"] = 2.0
_parameters["uuid_field"] = "uuid"
_parameters["height_field"] = ""

_parameters["statistical_outlier_remover"] = True
_parameters["roof_outlier_neighbors"] = 5
_parameters["roof_outlier_margin"] = 1.5

_parameters["ransac_outlier_remover"] = True
_parameters["ransac_outlier_margin"] = 3.0
_parameters["ransac_iterations"] = 250

_parameters["smoothing_max_iterations"] = 1000
_parameters["smoothing_relative_tolerance"] = 0.001

_parameters["naive_vegitation_filter"] = True
_parameters["data_directory"] = ""
_parameters["buildings_filename"] = "footprints.shp"
_parameters["pointcloud_directory"] = ""
_parameters["output_directory"] = ""


def default():
    return _parameters.copy()
