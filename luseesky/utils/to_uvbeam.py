from pyuvdata import uvbeam
from typing import Optional, NoReturn

def _flatten(arr, freq_size, th_size, ph_size):
    """
    Convert array with the shape (freq_size, th_size, ph_size) to a 2d array
    of shape (freq_size, th_size*ph_size) where theta increases faster than phi
    """
    return arr.reshape(freq_size, th_size*ph_size, order="F")

def _write_txt_power(
        beam: np.ndarray,  # must be in V
        frequencies: np.ndarray
        theta: np.ndarray,
        phi: np.ndarray,
        freq_axis: int = 0
        theta_axis: int = 1,
        phi_axis: int = 2,
        path: str = ""
) -> str:

    beam_2d = _flatten(beam, frequencies.size, theta.size, phi.size)
    for i, freq in enumerate(frequencies):
        np.savetxt(
            path+f"/{freq}.txt",
            np.column_stack((theta, phi, beam_2d[i]))
            header="Theta [rad] Phi [rad] Abs(V) [V] \n\n",
            comments="",
        )

    return path

def _write_txt_efield():
    raise (NotImplementedError)

def _delete_txt(path: str) -> NoReturn:
    # cleaning

def write_uvb(fname: str, beam: np.ndarray, frequencies: np.ndarray,
        beam_type: str = "power", outpath: str = ""): -> Optional[uvbeam.UVBeam]
    txtpath = _write_txt()
    uvb = uvbeam.UVBeam()
    uvb.read_cst_beam(
            filename=txtpath,
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
            x_orientation="north",
            reference_impedance=50)
    _delete_txt(txtpath)
    if len(outpath):
        uvb.write_beamfits(outpath, clobber=True)
        return None
    else:
        return uvb
