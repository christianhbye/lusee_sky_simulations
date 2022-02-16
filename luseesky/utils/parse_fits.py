from dataclasses import dataclass, field
from astropy.io import fits  # type: ignore
import matplotlib.pyplot as plt
import numpy as np
from typing import Any
import warnings

# import to_uvbeam which makes a tmp txt file and calls read cst beam


def mk_linspace(low: float, high: float, step: Any = 1) -> np.ndarray:
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
    num = int((high - low) // step + 1)
    return np.linspace(low, high, num=num)

def Efield_to_power(efield: np.ndarray, axis: int = 3) -> np.ndarray:
    """
    axis: the axis representing the x,yz-components. To be summed over.
    """
    return np.sum(np.abs(efield)**2, axis=axis)


@dataclass(frozen=True)
class Beam:
    __slots__ = ["fname"]
    fname: str
    E_field: np.ndarray = field(init=False)
    power: np.ndarray = field(init=False)
    frequencies: np.ndarray = field(init=False)
    theta: np.ndarray = field(init=False)
    phi: np.ndarray = field(init=False)

    def __post_init__(self):
        simfits = fits.open(self.fname)
        header = simfits[0].header
        self.E_field = fits[0].data + 1j * fits[2].data
        self.power = Efield_to_power(self.E_field, axis=3)
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
        

    def plot_power(self, freq: float):
        freq_idx = np.argmin(np.abs(self.frequencies-freq))
        plt.figure()
        plt.imshow(self.power[freq_idx], interpolation='none')  # set extent
        plt.colorbar()  # units of V squared
        plt.title(f"Power at $\\nu={self.frequencies[freq_idx]}$ MHz")
        plt.xlabel("$\\theta$ [deg]")
        plt.ylabel("$\\phi$ [deg]")
        plt.show()

    def plot_beamcuts(self, phi: float = 0): 
        phi_idx = np.argmin(np.abs(self.phi-phi))
        plt.figure()
        plt.imshow(self.power[:, :, phi_idx], interpolation='none')  # set extent
        plt.colorbar()  # units of V squared
        plt.title(f"Power at $\\phi={self.phi[phi_idx]}$ MHz")
        plt.xlabel("$\\nu$ [MHz]")
        plt.ylabel("$\\theta$ [deg]")
        plt.show()

    
    def to_uvbeam(self, save=True):
        uvb = mpkg.to_uvbeam(self.data)
        if save:
            uvb.write_beamfits()
        # import to_uvbeam which makes a tmp txt file and calls read cst beam

