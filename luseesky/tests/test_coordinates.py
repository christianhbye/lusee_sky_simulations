from luseesky.utils import coordinates as lcoords
import numpy as np

def test_sph2cart():
    E = np.array([1, 0, 0])  # radial vector along z
    th, ph = 0, 0
    rot = lcoords.sph2cart(th, ph)
    Ex, Ey, Ez = rot @ E
    assert np.allclose(rot.T, np.linalg.inv(rot))
    assert np.isclose(Ex, 0)
    assert np.isclose(Ey, 0)
    assert np.isclose(Ez, 1)
    assert np.isclose(np.linalg.norm(E), np.linalg.norm([Ex, Ey, Ez]))

    # x-direction
    E = np.array([0, 0, 1])
    th, ph = np.pi/2, -np.pi/2
    rot = lcoords.sph2cart(th, ph)
    Ex, Ey, Ez = rot @ E
    assert np.allclose(rot.T, np.linalg.inv(rot))
    assert np.isclose(Ex, 1)
    assert np.isclose(Ey, 0)
    assert np.isclose(Ez, 0)
    assert np.isclose(np.linalg.norm(E), np.linalg.norm([Ex, Ey, Ez]))

def test_cart2sph():
    r_hat = np.ones(3)/np.sqrt(3)  # radial unit vector @ (1, 1, 1)
    th_hat = np.array([0, 0, -1])  # polar unit vector @ (1, 0, 0)
    ph_hat = np.array([0, 1, 0])  # azimuthal unit vector @ (1, 0, 0)

    # radial unit vector from (1, 1, 1)
    E = r_hat
    th, ph = np.arccos(1/np.sqrt(3)), np.pi/4
    rot = lcoords.cart2sph(th, ph)
    Er, Eth, Eph = rot @ E
    assert np.allclose(rot.T, np.linalg.inv(rot))
    assert np.isclose(Er, 1)
    assert np.isclose(Eth, 0)
    assert np.isclose(Eph, 0)
    assert np.isclose(np.linalg.norm(E), np.linalg.norm([Er, Eth, Eph]))

    # polar unit vector at x=1, y=0, z=0
    E = th_hat
    th, ph = np.pi/2, 0
    rot = lcoords.cart2sph(th, ph)
    Er, Eth, Eph = rot @ E
    assert np.allclose(rot.T, np.linalg.inv(rot))
    assert np.isclose(Er, 0)
    assert np.isclose(Eth, 1)
    assert np.isclose(Eph, 0)
    assert np.isclose(np.linalg.norm(E), np.linalg.norm([Er, Eth, Eph]))

    # azimuthal unit vector at x=1, y=0, z=0
    E = ph_hat
    th, ph = np.pi/2, 0
    rot = lcoords.cart2sph(th, ph)
    Er, Eth, Eph = rot @ E
    assert np.allclose(rot.T, np.linalg.inv(rot))
    assert np.isclose(Er, 0)
    assert np.isclose(Eth, 0)
    assert np.isclose(Eph, 1)
    assert np.isclose(np.linalg.norm(E), np.linalg.norm([Er, Eth, Eph]))

