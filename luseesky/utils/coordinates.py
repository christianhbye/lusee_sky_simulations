import numpy as np


def sph2cart(theta: float, phi: float) -> np.ndarray:
    """
    Get the rotation matrix for transforming spherical coordinates
    to cartesian. Note that the matrix is orthogonal so the inverse transform
    is given by its transpose.
    """
    sin_th = np.sin(theta)
    cos_th = np.cos(theta)
    sin_ph = np.sin(phi)
    cos_ph = np.cos(phi)
    rot_matrix = np.array(
        [
            [sin_th * cos_ph, cos_th * cos_ph, -sin_ph],
            [sin_th * sin_ph, cos_th * sin_ph, cos_ph],
            [cos_th, -sin_th, 0],
        ]
    )
    return rot_matrix


def cart2sph(theta: float, phi: float) -> np.ndarray:
    """
    Get the rotation matrix for transforming cartesian coordinates to
    sphericals given sky angles. This is simply the transpose of the matrix that
    transforms sphericals to cartesian.
    """
    rot_matrix = sph2cart(theta, phi).T
    return rot_matrix
