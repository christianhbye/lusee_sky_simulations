from csv_to_txt import Beam
from misc import get_freq
from argparse import ArgumentParser
import numpy as np
from pathlib import Path
import pyuvdata

parser = ArgumentParser()
parser.add_argument("beam_dir", metavar="-b", type=str)
parser.add_argument("txt_out_dir", metavar="-t", type=str)
parser.add_argument("uvbeam_path", metavar="-uv", type=str)
args = parser.parse_args()

beam_path = Path(args.beam_dir)

beam_inphase = [p for p in beam_path.glob("*_in_phase.csv") if "perpendicular"
        not in p.name]
beam_outphase = [p for p in beam_path.glob("*_out_of_phase.csv") if "perpendicular" not in p.name]

out_paths = []
frequencies = []
for fname in beam_inphase:
    freq = get_freq(fname.name, unit="MHz")
    frequencies.append(freq)
    out_path = Path(args.txt_out_dir) / Path(f"Freq{freq}MHz_in_phase.txt")
    out_paths.append(str(out_path))
    beam = Beam(str(fname))
    beam.to_txt(out_path)

out_paths = sorted(out_paths, key=lambda x: frequencies[out_paths.index(x)])
frequencies = sorted(frequencies)

uvb = pyuvdata.uvbeam.UVBeam()
uvb.read_cst_beam(
       filename=out_paths,
       beam_type="power",
       feed_pol="x",
       rotate_pol=False,
       frequency=frequencies,
       telescope_name="lusee-night",
       feed_name="lusee",
       feed_version="1.0",
       model_name="monopole-in-phase",
       model_version="1.0",
       history="002KMR",
       x_orientation="north",
       reference_impedance=50
       )

uvb.write_beamfits(args.uvbeam_path, clobber=True)
