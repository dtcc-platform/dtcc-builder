import _pybuilder


def generateElevationModel(point_cloud, resolution, included_classifications=[]):
    included_classifications = list(map(int, included_classifications))
    grid_field = _pybuilder.GenerateElevationModel(
        point_cloud, resolution, included_classifications
    )
    return grid_field


def smoothElevationModel(grid_field, num_passes):
    smoothed_field = _pybuilder.SmoothElevation(grid_field, num_passes)
    return smoothed_field

def getBounds(grid_field):
    p = grid_field.Grid.BoundingBox.P
    q = grid_field.Grid.BoundingBox.Q
    return (p.x,p.y,q.x,q.y)
