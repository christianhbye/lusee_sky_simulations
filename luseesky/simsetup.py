from pathlib import Path
import yaml

def gen_uvbdict(ant_path: str, outpath: str = "./sim_files/uvbeam.yaml"):
    """
    Generate the uvbeam.yaml file.
    """
    uvbdict = {
            "beam_paths": {
                0: str(Path(ant_path).resolve()),
                1: {
                    "type": "airy",
                    "diameter": 16
                    }
                },
            "freq_interp_kind": "linear",
            "telescope_location": "(0., 180., 0.)",
            "world": "moon",
            "telescope_name": "LuSEE-Night",
            "x_orientation": "north"
            }
    with open(outpath, "w") as f:
        yaml.safe_dump(uvbdict, f)

def gen_obsparams(
        ant_model: str,
        nsrcs: int = 500,
        outpath: str = "./sim_files/obsparam.yaml"
):
    """
    Generate the obsparam dict.
    """
    obsparams = {
            "filing": {
                "outdir": "./results/results",
                "outfile_name": f"{ant_model}_uvsim",
                "output_format": "uvfits"
                },
            "freq": {
                "Nfreqs": 50,
                "start_freq": 1000000.0,
                "end_freq": 50000000.0,
                "channel_width": 1000000.0
                },
            "sources": {
                "catalog": f"./skymodels/pygds16_{nsrcs}srcs.txt"
                },
            "telescope": {
                "array_layout": "./sim_files/layout.csv",
                "telescope_config_name": "./sim_files/uvbeam.yaml"
                },
            "time": {
                "Ntimes": 7680,
                "start_time": 2459630.,
                "duration_days": 30
                },
            "select": {
                "bls": "[(0, 0)]"
                }
            }
    with open(outpath, "w") as f:
        yaml.safe_dump(obsparams, f)

