import _pybuilder


class ElevationModel:
    def __init__(self, point_cloud, resolution, included_classifications=[]):
        self.included_classifications = list(map(int, included_classifications))
        self.resolution = resolution
        self._grid_field = _pybuilder.GenerateElevationModel(
            point_cloud, resolution, included_classifications
        )
        p = self._grid_field.Grid.BoundingBox.P
        q = self._grid_field.Grid.BoundingBox.Q
        self.bounds = (p.x, p.y, q.x, q.y)

    def smooth_elevation_model(self, num_passes):
        self._grid_field = _pybuilder.SmoothElevation(self._grid_field, num_passes)
