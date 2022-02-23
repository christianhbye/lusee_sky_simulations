from argparse import ArgumentParser
from luseesky.utils.parse_fits import Beam
from pathlib import Path

parser = ArgumentParser()
parser.add_argument("dir", type=str)
args = parser.parse_args()

DIR = args.dir

UVDIR = "../uvbeams/" + DIR[len("../../"):]
assert "uvbeams/Ant" in UVDIR  # checking that joining the dirname works
uvp = Path(UVDIR)
if not uvp.is_dir():
    try:
        uvp.mkdir(parents=True)
    except (FileExistsError):
        print("The path exists but it is not a dir ...?")

for beamfile in Path(DIR).iterdir():
    if not beamfile.suffix == ".fits":
        continue
    UVNAME = beamfile.name[:-(len(".fits"))] + ".uvbeam"
    beam = Beam(str(beamfile))
    uvb = beam.to_uvbeam()
    try:
        uvb.write_beamfits(UVDIR+"/"+UVNAME, clobber=False)
    except(OSError) as e:
        print(e)
        continue
