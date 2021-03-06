# from luseesky import __path__
import numpy as np
from pathlib import Path
from pyuvsim import uvsim  # type: ignore
from luseesky.simsetup import gen_uvbdict, gen_obsparams

# LPATH = __path__[0]
LPATH = (
    "/home/christian/Documents/research/lusee/lusee_sky_simulations/"
    "luseesky"
)


def run(outpath: str):
    uvd = uvsim.run_uvsim(
        LPATH + "/sim_files/obsparam.yaml", return_uv=True, quiet=False
    )
#    uvd.data_array = uvd.data_array.real

    save_dict = {
        "data": uvd.data_array,
        "freq": uvd.freq_array,
        "lst": uvd.lst_array,
    }

    np.savez(outpath, **save_dict)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("dir", type=str)  # path to uvbeam of antenna
    args = parser.parse_args()

    DIR = args.dir

    RESDIR = LPATH + "/results/" + DIR[len("uvbeams/AntennaSimResults/") :]
    assert "results/004" in RESDIR  # checking that joining the dirname works
    rp = Path(RESDIR)
    if not rp.is_dir():
        try:
            rp.mkdir(parents=True)
        except (FileExistsError):
            print("The path exists but it is not a dir ...?")

    for beamfile in Path(DIR).iterdir():
        if not beamfile.suffix == ".uvbeam":
            continue
        RESNAME = beamfile.name[: -len(".uvbeam")]
        gen_uvbdict(str(beamfile))
        gen_obsparams(str(beamfile))
        run(RESDIR + "/" + RESNAME)
