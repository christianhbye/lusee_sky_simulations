import csv_to_txt
from argparse import ArgumentParser
import numpy as np
from pathlib import Path
import pyuvdata

parser = ArgumentParser()
parser.add_argument("beam_dir", type=str)
args = parser.parse_args()

beam_path = Path(args.beam_dir)

frequencies = np.unique(
    [
        float(fname.name.split("Freq")[-1].split("MHz")[0])
        for fname in beam_path.iterdir()
    ]
)

beam_inphase = beam_path.glob("*_in_phase.csv")
beam_outphase = beam_path.glob("*_out_of_phase.csv")
