from dataclasses import dataclass, field
from astropy.io import fits
import numpy as np
import warnings


def mk_linspace(low: float, high: float, step: float = 1) -> np.ndarray:
    """
    Make a linspace given low, high, step. This avoids the stability
    issues in np.arange(low, high, step) when low >> step (see numpy doc).
    If high is not a multiple of steps away from low, a warning will be raised.
    """
    if not np.isclose((high - low) // step, (high - low) / step, atol=1e-4):
        warnings.warn(
            "'high' is not a multiple of 'step' away from 'low',\
                'step' will be changed.",
            UserWarning,
        )
    num = (high - low) // step + 1
    return np.linspace(low, high, num=num)


@dataclass(frozen=True)
class Beam:
    __slots__ = ["fname", "frequencies", ""]
    fname: str
    E_field: np.ndarray = field(init=False)
    frequencies: np.ndarray = field(init=False)
    theta: np.ndarray = field(init=False)
    phi: np.ndarray = field(init=False)

    def __post_init__(self):
        simfits = fits.open(self.fname)
        header = simfits[0].header
        self.E_field = fits[0].data + 1j * fits[2].data
        self.frequencies = mk_linspace(
            header["freq_start"], header["freq_end"], step=header["freq_step"]
        )
        self.theta = mk_linspace(
            header["theta_start"],
            header["theta_end"],
            step=header["theta_step"],
        )
        self.phi = mk_linspace(
            header["phi_start"], header["phi_end"], step=header["phi_step"]
        )
