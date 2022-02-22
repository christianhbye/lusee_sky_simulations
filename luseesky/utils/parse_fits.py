from astropy.io import fits  # type: ignore
from dataclasses import dataclass, field
import matplotlib.pyplot as plt  # type: ignore
import numpy as np
from pathlib import Path
from pyuvdata import uvbeam  # type: ignore
from typing import Any, Tuple
import warnings


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

    returns beam in V
    """

    return np.sqrt(np.sum(np.abs(efield) ** 2, axis=axis))


@dataclass
class Beam:
    fname: str
    E_field: np.ndarray = field(init=False)
    power: np.ndarray = field(init=False)
    frequencies: np.ndarray = field(init=False)
    theta: np.ndarray = field(init=False)
    phi: np.ndarray = field(init=False)

    def __post_init__(self):
        simfits = fits.open(self.fname)
        header = simfits[0].header
        self.E_field = simfits[0].data + 1j * simfits[1].data
        simfits.close()
        self.E_field /= 1e3  # convert mV to V
        self.frequencies = mk_linspace(
            header["freq_start"], header["freq_end"], step=header["freq_step"]
        )  # in MHz
        self.theta = mk_linspace(
            header["theta_start"],
            header["theta_end"],
            step=header["theta_step"],
        )
        self.phi = mk_linspace(
            header["phi_start"], header["phi_end"], step=header["phi_step"]
        )
        if np.allclose(self.E_field[:, :, 0], self.E_field[:, :, -1]):
            assert np.isclose(self.phi[-1] - self.phi[0], 360)
            self.E_field = self.E_field[:, :, :-1]  # drop phi = 360 deg
            self.phi = self.phi[:-1]
        self.power = Efield_to_power(self.E_field, axis=3)

    def plot_power(self, freq: float):
        freq_idx = np.argmin(np.abs(self.frequencies - freq))
        plt.figure()
        plt.imshow(
            self.power[freq_idx],
            interpolation="none",
            aspect="auto",
            extent=[
                self.phi.min(),
                self.phi.max(),
                self.theta.max(),
                self.theta.min(),
            ],
        )
        plt.colorbar(label="Power [V]")
        plt.title(
            "Power at $\\nu={:.0f}$ MHz".format(self.frequencies[freq_idx])
        )
        plt.xlabel("$\\phi$ [deg]")
        plt.ylabel("$\\theta$ [deg]")
        plt.show()

    def plot_beamcuts(self, phi: float = 0):
        phi_idx = np.argmin(np.abs(self.phi - phi))
        plt.figure()
        plt.imshow(
            self.power[:, :, phi_idx],
            interpolation="none",
            aspect="auto",
            extent=[
                self.theta.min(),
                self.theta.max(),
                self.frequencies.max(),
                self.frequencies.min(),
            ],
        )
        plt.title("Power at $\\phi={:.0f}$ deg".format(self.phi[phi_idx]))
        plt.ylabel("$\\nu$ [MHz]")
        plt.xlabel("$\\theta$ [deg]")
        plt.colorbar(label="Power [V]")
        plt.show()

    def _flatten(
        self, beam_type: str = "power"
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Convert array with the shape (freq_size, th_size, ph_size) to a
        2d-array of shape (freq_size, th_size*ph_size) where theta increases
        faster than phi
        """
        if beam_type == "power":
            arr = np.copy(self.power)
        elif beam_type == "E_field":
            raise NotImplementedError("Might add E_field support later.")
        else:
            raise ValueError("beam_type must be 'power' or 'E_field'")
        flat_beam = arr.reshape(
            self.frequencies.size, self.theta.size * self.phi.size, order="F"
        )
        flat_theta = np.tile(self.theta, self.phi.size)
        flat_phi = (
            np.tile(self.phi, self.theta.size)
            .reshape(self.phi.size, self.theta.size, order="F")
            .flatten(order="C")
        )
        return flat_beam, flat_theta, flat_phi

    def _write_txt_power(self, path: str = ".", verbose: bool = False) -> str:
        beam2d, th2d, ph2d = self._flatten()
        savepath = path + "/tmp"
        Path(savepath).mkdir()
        if verbose:
            print(f"Saving {len(beam2d)} files to {savepath}")
        for i, freq in enumerate(self.frequencies):
            np.savetxt(
                savepath + f"/{freq}.txt",
                np.column_stack((th2d, ph2d, beam2d[i])),
                header="Theta [deg] Phi [deg] Abs(V) [V] \n\n",
                comments="",
            )
        return savepath

    @staticmethod
    def _delete_txt(path: str, verbose: bool = False):
        for f in Path(path).iterdir():
            assert f.suffix == ".txt"
            f.unlink()  # delete file
        if verbose:
            print("Deleting files.")
        Path(path).rmdir()
        if verbose:
            print(f"Remove directory {path}.")

    def to_uvbeam(
        self, beam_type: str = "E_field", verbose: bool = False
    ) -> uvbeam.UVBeam:
        uvb = uvbeam.UVBeam()
        if beam_type == "power":
            if verbose:
                print("Making UVBeam object from power beam.")
            txtpath = self._write_txt_power(verbose=verbose)
            txtfiles = [str(child) for child in Path(txtpath).iterdir()]
            frequencies = [
                1e6 * float(Path(f).name[: -len(".txt")]) for f in txtfiles
            ]
            txtfiles = sorted(
                txtfiles, key=lambda x: frequencies[txtfiles.index(x)]
            )
            frequencies = sorted(frequencies)
            uvb.read_cst_beam(
                filename=txtfiles,
                beam_type="power",
                feed_pol="x",
                rotate_pol=False,
                frequency=frequencies,
                telescope_name="lusee-night",
                feed_name="lusee",
                feed_version="1.0",
                model_name="monopole",
                model_version="1.0",
                history="003",
                reference_impedance=50,
            )
            uvb.interpolation_function = "az_za_simple"
            self._delete_txt(txtpath, verbose=verbose)
        elif beam_type == "E_field":
            if verbose:
                print("Making UVBeam object from E-field beam.")
            uvb.filename = [Path(self.fname).name]
            uvb._filename.form = (1,)
            uvb.telescope_name = "lusee-night"
            uvb.feed_name = "lusee"
            uvb.feed_version = "1.0"
            uvb.model_name = "monopole"
            uvb.model_version = "1.0"
            uvb.history = "003" + uvb.pyuvdata_version_str
            uvb.reference_impedance = 50.0
            uvb.Naxes_vec = 2
            uvb.Ncomponents_vec = 2
            uvb.feed_array = np.array(["x", "y"])
            uvb.Nfeeds = uvb.feed_array.size
            uvb._set_efield()
            uvb.data_normalization = "physical"
            uvb.antenna_type = "simple"
            uvb.Nfreqs = self.frequencies.size
            uvb.Nspws = 1
            uvb.freq_array = self.frequencies.reshape(1, -1) * 1e6
            uvb.bandpass_array = np.zeros_like(uvb.freq_array)
            uvb.spw_array = np.array([0])
            uvb.pixel_coordinate_system = "az_za"
            uvb._set_cs_params()
            uvb.axis1_array = np.radians(self.phi)
            uvb.Naxes1 = uvb.axis1_array.size
            uvb.axis2_array = np.radians(self.theta)
            uvb.Naxes2 = uvb.axis2_array.size
            uvb.data_array = np.zeros(
                uvb._data_array.expected_shape(uvb), dtype="complex128"
            )
            uvb.basis_vector_array = np.zeros(
                (uvb.Naxes_vec, uvb.Ncomponents_vec, uvb.Naxes2, uvb.Naxes1)
            )
            uvb.basis_vector_array[0, 0] = 1.0
            uvb.basis_vector_array[1, 1] = 1.0
            # data_array: [x,y], 0, [feed/pol], freq, theta, phi
            uvb.data_array[0, 0, 0] = self.E_field[:, :, :, 0]
            uvb.data_array[1, 0, 1] = self.E_field[:, :, :, 1]
            uvb.bandpass_array[0] = 1
            uvb.check(check_extra=True, run_check_acceptability=False)
            uvb.interpolation_function = "az_za_simple"
        else:
            raise ValueError("beam_type must be 'power' or 'E_field'")
        return uvb
