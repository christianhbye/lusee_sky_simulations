from luseesky.utils import to_uvbeam as ltu
import numpy as np


def test_flatten():
    # generate a beam
    freq_size, th_size, ph_size = 50, 91, 180
    rng = np.random.default_rng(seed=42)
    data = rng.normal(size=(freq_size, th_size, ph_size))
    flat = ltu._flatten(data)  # make 2d
    for phi in [0, 5, 45, 137]:  # fixed phi and frequency (=30)
        assert np.allclose(
                data[30, :, phi],
                flat[30, phi*th_size:(phi+1)*th_size]
            )
    for theta in [13, 27, 50, 88]:  # fixed theta and frequency (=5)
        assert np.allclose(data[5, theta, :], flat[5, theta::th_size])
