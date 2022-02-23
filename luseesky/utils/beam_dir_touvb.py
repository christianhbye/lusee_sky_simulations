from luseesky.utils.parse_fits import Beam
from pathlib import Path

DIR = ("AntennaSimResults/003_Freq1-50MHz_Delta1MHz_AntennaLength6m_Antenna"
        "Angle30deg_LanderHeight2m")

UVDIR = "../uvbeams/" + DIR
uvp = Path(UVDIR)
if not uvp.is_dir():
    try:
        uvp.mdkir(parents=True)
    except (FileExistsError):
        print("The path exists but it is not a dir ...?")

for beamfile in Path(DIR).iterdir():
    if beamfile.is_dir():
        continue
    UVNAME = beamfile.name[:-(len(".fits"))] + ".uvbeam"
    beam = Beam(beamfile.name)
    uvb = beam.to_uvbeam()
    try:
        uvb.write_beamfits(UVDIR+UVNAME, clobber=False)
    except(OSError) as e:
        print(e)
        continue
