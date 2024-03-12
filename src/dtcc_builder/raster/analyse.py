import numpy as np
from scipy.signal import convolve2d
from scipy.ndimage import convolve
from dtcc_model import Raster
from dtcc_builder.register import register_model_method


@register_model_method
def slope_aspect(dem: Raster) -> tuple[Raster, Raster]:
    cell_size = dem.cell_size[0]
    kernel_x = np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]]) / (8.0 * cell_size)
    kernel_y = np.array([[1, 2, 1], [0, 0, 0], [-1, -2, -1]]) / (8.0 * cell_size)

    # dzdx = convolve(dem.data, kernel_x)
    dzdx = convolve2d(dem.data, kernel_x, mode="same", boundary="symm")

    # dzdy = convolve(dem.data, kernel_y)
    dzdy = convolve2d(dem.data, kernel_y, mode="same", boundary="symm")
    slope = np.arctan(np.sqrt(dzdx**2 + dzdy**2))
    aspect = np.arctan2(-dzdy, dzdx)

    slope_raster = dem.copy(no_data=True)
    slope_raster.data = slope

    aspect_raster = dem.copy(no_data=True)
    aspect_raster.data = aspect

    return slope_raster, aspect_raster


@register_model_method
def TRI(dem: Raster) -> Raster:
    """
    Compute the Terrain Roughness Index (TRI) of a DEM using optimized array operations.

    Parameters:
    - dem: a Raster object representing the DEM.

    Returns:
    - A Raster object representing the TRI.
    """

    window_size = 3
    # Create a kernel that has ones in all positions but the center.
    kernel = np.ones((window_size, window_size))

    dem2 = dem.data**2
    s = convolve(dem.data, kernel)
    t = convolve(dem2, kernel)

    squared_diff = t - 2 * s * dem.data + (window_size * window_size) * dem2

    # Compute the squared difference between each cell and its neighbors.
    # squared_diff = (
    #     convolve2d(dem.data**2, kernel, mode="same", boundary="symm")
    #     - 2 * convolve2d(dem.data, kernel, mode="same", boundary="symm")
    #     + (window_size * window_size) * dem.data**2
    # )

    # Compute the TRI.
    tri = np.sqrt(squared_diff)

    tri_raster = dem.copy(no_data=True)
    tri_raster.data = tri

    return tri_raster


@register_model_method
def TPI(dem: Raster, window_size=3):
    """Compute the Topographic Position Index (TPI) of a DEM."""

    if (window_size % 2) == 0 or window_size < 2:
        raise ValueError("Window size must be odd and greater than 1.")

    # Define averaging kernel
    kernel = np.ones((window_size, window_size)) / (window_size**2)

    # Compute the mean elevation of the neighborhood
    mean_elevation = convolve2d(dem.data, kernel, mode="same", boundary="symm")

    # Compute TPI
    tpi = dem.data - mean_elevation
    tpi_raster = dem.copy(no_data=True)
    tpi_raster.data = tpi

    return tpi_raster


@register_model_method
def VRM(dem: Raster, window_size=3):
    """Compute the Vector Ruggedness Measure (VRM) of a DEM."""

    if (window_size % 2) == 0 or window_size < 2:
        raise ValueError("Window size must be odd and greater than 1.")
    slope, aspect = slope_aspect(dem)
    slope = slope.data
    aspect = aspect.data
    # Convert slope and aspect to 3D unit vectors
    x = np.cos(slope) * np.cos(aspect)
    y = np.cos(slope) * np.sin(aspect)
    z = np.sin(slope)

    # Define averaging kernel
    kernel = np.ones((window_size, window_size))

    # Average the unit vectors over the window
    x_avg = convolve2d(x, kernel, mode="same", boundary="symm") / (window_size**2)
    y_avg = convolve2d(y, kernel, mode="same", boundary="symm") / (window_size**2)
    z_avg = convolve2d(z, kernel, mode="same", boundary="symm") / (window_size**2)

    # Compute VRM
    magnitude_avg = np.sqrt(x_avg**2 + y_avg**2 + z_avg**2)
    vrm = np.sqrt(1 - (magnitude_avg / window_size) ** 2)

    vrm_raster = dem.copy(no_data=True)
    vrm_raster.data = vrm

    return vrm_raster
