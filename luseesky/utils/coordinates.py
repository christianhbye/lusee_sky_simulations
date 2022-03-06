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


# def sph2cart_array(theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
#     """
#     Calling sph2cart for a range of thetas and phis.
#     """
#     master_rot = np.empty((theta.size, phi.size, 3, 3))
#     for i, th in enumerate(theta):
#         for j, ph in enumerate(phi):
#             master_rot[i, j] = sph2cart(th, ph)
#     return master_rot
#
# def cart2sph_array(theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
#     """
#     Calling cart2sph for a range of thetas and phis.
#     """
#     master_rot = np.empty((theta.size, phi.size, 3, 3))
#     for i, th in enumerate(theta):
#         for j, ph in enumerate(phi):
#             master_rot[i, j] = cart2sph(th, ph)
#     return master_rot
