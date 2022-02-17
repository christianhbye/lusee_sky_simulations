from luseesky.utils import parse_fits as lpf
import numpy as np
from pathlib import Path
import pytest


def test_mk_linspace():
    lo, hi, step = -1, 5.5, 0.5
    arr = lpf.mk_linspace(lo, hi, step=step)
    assert np.isclose(lo, arr.min())
    assert np.isclose(hi, arr.max())
    assert np.allclose(np.diff(arr), step*np.ones_like(np.diff(arr)))

    # should raise a warning
    lo, hi, step = 1, 3.2, 0.5
    with pytest.warns(UserWarning):
        arr = lpf.mk_linspace(lo, hi, step)
        # should still be true:
        assert np.isclose(lo, arr.min())
        assert np.isclose(hi, arr.max())

def test_flatten():
    assert Path("AntennaSimResults").exists()
    test_beam = lpf.Beam(
            "AntennaSimResults/003_Freq1-50MHz_Delta1MHz_AntennaLength6m"
            "_AntennaAngle30deg_LanderHeight2m/RadiatedElectricField_Antenna"
            "Length6m_AntennaAngle30deg_LanderHeight2m_monopole.fits"
            )
    flat_beam = test_beam._flatten(beam_type="power")
    th_size = test_beam.power.shape[1]
    ph_size = test_beam.power.shape[2]
    for phi in [0, 5, 45, 137]:  # fixed phi and frequency (=30)
        assert np.allclose(
                test_beam.power[30, :, phi],
                flat_beam[30, phi*th_size:(phi+1)*th_size]
            )
    for theta in [13, 27, 50, 88]:  # fixed theta and frequency (=5)
        assert np.allclose(
                test_beam.power[5, theta, :],
                flat_beam[5, theta::th_size]
            )

