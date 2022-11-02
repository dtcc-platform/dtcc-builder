import _pybuilder

def generateElevationModel(point_cloud,resolution, included_classifications=[]):
    included_classifications = list(map(int,included_classifications))
    grid_field = _pybuilder.GenerateElevationModel(point_cloud,resolution,included_classifications)
    return grid_field

def smoothElevationModel(grid_field,num_passes):
    smoothed_field = _pybuilder.SmoothElevation(grid_field,num_passes)
    return smoothed_field