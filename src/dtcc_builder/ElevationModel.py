from google.protobuf.json_format import MessageToJson

from dtcc_builder import _pybuilder


class ElevationModel:
    def __init__(self, point_cloud=None, resolution=1, included_classifications=[]):
        self.included_classifications = list(map(int, included_classifications))
        self.resolution = resolution
        self.bounds = (0, 0, 0, 0)
        if point_cloud is not None:
            self._grid_field = _pybuilder.GenerateElevationModel(
                point_cloud._builder_pc, resolution, included_classifications
            )
            self.calc_bounds()

    def calc_bounds(self):
        p = self._grid_field.Grid.BoundingBox.P
        q = self._grid_field.Grid.BoundingBox.Q
        self.bounds = (p.x, p.y, q.x, q.y)

    def smooth_elevation_model(self, num_passes=5):
        if num_passes < 1:
            return
        self._grid_field = _pybuilder.SmoothElevation(self._grid_field, num_passes)

    def mean(self):
        return _pybuilder.MeanElevation(self._grid_field)

    def min(self):
        return _pybuilder.MinElevation(self._grid_field)

    def max(self):
        return _pybuilder.MaxElevation(self._grid_field)

    # TODO: Implement this
    def to_protobuf(self):
        return ""
