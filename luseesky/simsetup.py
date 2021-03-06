# from luseesky import __path__
from pathlib import Path
import yaml  # type: ignore

# LPATH = __path__[0]
LPATH = (
    "/home/christian/Documents/research/lusee/lusee_sky_simulations/"
    "luseesky"
)


def gen_uvbdict(
    ant_path: str,
    telescope_coords: str = "(0., 180., 0.)",
    outpath: str = LPATH + "/sim_files/uvbeam.yaml"
):
    """
    Generate the uvbeam.yaml file.

    Parameters
        ant_path: str, path to the .uvbeam file of the antenna
        telescope_coords: str, moon coordinates (lon, lat, alt) of telescope
        outpath: str, path + filename of the uvbeam.yaml file that goes into
        pyuvsim
    """
    uvbdict = {
        "beam_paths": {
            0: str(Path(ant_path).resolve()),
            1: {"type": "airy", "diameter": 16},
        },
        "freq_interp_kind": "linear",
        "telescope_location": f"{telescope_coords}",
        "world": "moon",
        "telescope_name": "LuSEE-Night",
        "x_orientation": "north",
    }
    with open(outpath, "w") as f:
        yaml.safe_dump(uvbdict, f)


def gen_obsparams(
    ant_model: str,
    nside: int = 64,
    outpath: str = LPATH + "/sim_files/obsparam.yaml",
):
    """
    Generate the obsparam dict.

    Parameters:
        ant_model: str, some defining characteristic of antenna. This defines
        the outfile name so picking something unique lessens the chances of
        overwriting another result!
        nside: int, nside of healpix map of sky
        outpath: str, path + filename of obsparam.yaml file that gets fed into
        pyuvsim.
    """
    obsparams = {
        "filing": {
            "outdir": f"{LPATH}/results/results",
            "outfile_name": f"{ant_model}_uvsim",
            "output_format": "uvfits",
        },
        "freq": {
            "Nfreqs": 50,
            "start_freq": 1000000.0,
            "end_freq": 50000000.0,
            "channel_width": 1000000.0,
        },
        "sources": {
            "catalog": f"{LPATH}/skymodels/pygdsm16_nside{nside}.txt"
            },
        "telescope": {
            "array_layout": f"{LPATH}/sim_files/layout.csv",
            "telescope_config_name": f"{LPATH}/sim_files/uvbeam.yaml",
        },
        "time": {"Ntimes": 256, "start_time": 2459630.0, "duration_days": 28},
        "select": {"bls": "[(0, 0)]"},
    }
    with open(outpath, "w") as f:
        yaml.safe_dump(obsparams, f)
