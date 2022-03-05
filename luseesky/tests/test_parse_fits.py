from luseesky.utils import parse_fits as lpf
import numpy as np
import pytest


def test_mk_linspace():
    lo, hi, step = -1, 5.5, 0.5
    arr = lpf.mk_linspace(lo, hi, step=step)
    assert np.isclose(lo, arr.min())
    assert np.isclose(hi, arr.max())
    assert np.allclose(np.diff(arr), step * np.ones_like(np.diff(arr)))

    # should raise a warning
    lo, hi, step = 1, 3.2, 0.5
    with pytest.warns(UserWarning):
        arr = lpf.mk_linspace(lo, hi, step)
        # should still be true:
        assert np.isclose(lo, arr.min())
        assert np.isclose(hi, arr.max())


def test_flatten():
    test_beam = lpf.Beam("luseesky/tests/beam.fits")
    flat_beam, flat_th, flat_ph = test_beam._flatten(beam_type="power")
    th_size = test_beam.power.shape[1]
    ph_size = test_beam.power.shape[2]
    for phi in [0, 5, 45, 137]:  # fixed phi and frequency (=30)
        assert np.allclose(
            test_beam.power[30, :, phi],
            flat_beam[30, phi * th_size : (phi + 1) * th_size],
        )
    for theta in [13, 27, 50, 88]:  # fixed theta and frequency (=5)
        assert np.allclose(
            test_beam.power[5, theta, :], flat_beam[5, theta::th_size]
        )
    assert all(
        [
            np.allclose(
                flat_th[:th_size], flat_th[i * th_size : (i + 1) * th_size]
            )
            for i in range(ph_size)
        ]
    )  # flat theta repeats it self
    assert all(
        [
            np.allclose(flat_ph[::th_size], flat_ph[i::th_size])
            for i in range(th_size)
        ]
    )

def test_sph2cart():
    E = np.array([1, 0, 0])  # radial vector along z
    Ex, Ey, Ez = sph2cart(*E) @ E
    assert np.allclose(sph2cart(*E).T, np.linalg.inv(sph2cart(*E))
    assert np.isclose(Ex, 0)
    assert np.isclose(Ey, 0)
    assert np.isclose(Ez, 1)
    assert np.isclose(np.linalg.norm(E), np.linalg.norm([Ex, Ey, Ez])

    E = np.array([1, 1, 1])
    Ex, Ey, Ez = sph2cart(*E) @ E
    assert np.allclose(sph2cart(*E).T, np.linalg.inv(sph2cart(*E))
    assert np.isclose(Ex, 1 - np.sqrt(2)/2)
    assert np.isclose(Ey, 1 + np.sqrt(2)/2)
    assert np.isclose(Ez, 0)
    assert np.isclose(np.linalg.norm(E), np.linalg.norm([Ex, Ey, Ez])
    

def test_cart2sph():
    E = np.array([1, 1, 1])
    Er, Eth, Eph = cart2sph(*E) @ E
    assert np.allclose(cart2sph(*E).T, np.linalg.inv(cart2sph(*E))
    assert np.isclose(Er, 1 + np.sqrt(2)/2)
    assert np.isclose(Eth, 1 - np.sqrt(2)/2)
    assert np.isclose(Eph, 0)
    assert np.isclose(np.linalg.norm(E), np.linalg.norm([Er, Eth, Eph])

    E = [1, -1, 0]
    Er, Eth, Eph = cart2sph(*E) @ E
    assert np.allclose(cart2sph(*E).T, np.linalg.inv(cart2sph(*E))
    assert np.isclose(Er, np.sqrt(2))
    assert np.isclose(Eth, 0)
    assert np.isclose(Eph, 0)
    assert np.isclose(np.linalg.norm(E), np.linalg.norm([Er, Eth, Eph])

def test_to_sphericals_cartesian():
    """
    Check that both the real and imaginary parts are transformed correctly
    for all frequencies and sky coordinates
    """
    test_beam = lpf.Beam("luseesky/tests/beam.fits")
    assert test_beam.beam_coords == "cartesian"
    Ex, Ey, Ez = test_beam.E_field
    test_beam.to_sphericals()
    assert test_beam.beam_coords == "sphericals"
    E_sph = test_beam.E_field
    assert np.allclose(E_sph.real, cart2sph(Ex.real, Ey.real, Ez.real)
    assert np.allclose(E_sph.imag, cart2sph(Ex.imag, Ey.imag, Ez.imag)
    # convert back
    test_beam.to_cartesian()
    assert test_beam.beam_coords == "cartesian"
    assert np.allclose(np.array([Ex, Ey, Ez]), test_beam.E_field)

