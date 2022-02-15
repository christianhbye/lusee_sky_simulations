from luseesky.utils import parse_fits as lpf
import numpy as np
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
