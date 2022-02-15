from dataclasses import dataclass, field
import fitsio
import numpy as np
import warnings

def mk_linspace(
        low: int | float,
        high: int | float,
        step: int | float = 1
    ) -> np.ndarray:
    """
    Make a linspace given low, high, step. This avoids the stability
    issues in np.arange(low, high, step) when low >> step (see numpy doc).
    If high is not a multiple of steps away from low, a warning will be raised.
    """
    if not np.isclose((high-low)//step, (high-low)/step, atol=1e-4):
        warnings.warn(
                "'high' is not a multiple of 'step' away from 'low', 'step'
                will be changed.",
                UserWarning
            )
    num = (high-low)//step + 1
    return np.linspace(low, high, step)

@dataclass(frozen=True)
class Beam:
    __slots__ = ["fname", "frequencies", ""]
    fname: str
    E_field: np.ndarray = field(init=False)
    frequencies: np.ndarray = field(init=False)
    theta: np.ndarray = field(init=False)
    phi: np.ndarray = field(init=False)

    def __post_init__(self):
        header = fitsio.read_header(self.fname)
        fits = fitsio.FITS(self.fname, "r")
        self.E_field = fits[0].read() + 1j*fits[1].read()
        self.frequencies = mk_linspace(
                header["freq_start"],
                header["freq_end"],
                step=header["freq_step"]
            )    
        self.theta = mk_linspace(
                header["theta_start"],
                header["theta_end"],
                step=header["theta_step"]
            )
        self.phi = mk_linspace(
                header["phi_start"],
                header["phi_end"],
                step=header["phi_step"]
            )

