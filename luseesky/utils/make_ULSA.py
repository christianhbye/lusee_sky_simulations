import h5py
import numpy as np
from ULSA.sky_map.produce_absorbed_sky_map import absorption_JRZ

FREQS = 25
NSIDE = 32
INDEX_TYPE = "constant_index"

# default params
DISTANCE = 50  # kpc
PAR_408 = np.array([43.09932087,3.40820378,0.46230938,1.11894379,1.2344227])
RAW_DIFFUSE = False
FDIRS = None
DEF_PARS = True
IN_SPEC_IDX = None
CRIT_DI = False
OUT_FREE = False

sky_map = absorption_JRZ(
        v=freq,
        nside=NSIDE,
        index_type=INDEX_TYPE,
        distance=DISTANCE,
        using_raw_diffuse=RAW_DIFFUSE,
        v_file_dir=FDIRS,
        using_default_params=DEF_PARS,
        input_spectral_index=IN_SPEC_IDX,
        params_408=PAR_408,
        critical_dis=CRIT_DI,
        output_absorp_free_skymap=OUT_FREE
)

with h5py.File("../skymodels/skymap.h5", "w") as hf:
    hf.create_dataset("freqs_MHz", data=FREQS)
    hf.create_datset("nside", data=NSIDE)
    hf.create_dataset("index_type", data=INDEX_TYPE)
    hf.create_dataset("distance_kpc", data=DISTANCE)
    hf.create_dataset("params_408", data=PAR_408)
    hf.create_dataset("raw_diffuse", data=RAW_DIFFUSE)
    hf.create_dataset("freq_file_dirs", data=FDIRS)
    hf.create_dataset("default_params", data=DEF_PARS)
    hf.create_dataset("input_spectral_idx", data=IN_SPEC_IDX)
    hf.create_dataset("critical_dist", data=CRIT_DI)
    hf.create_dataset("output_abs_free", data=OUT_FREE)
    hf.create_dataset("skymap", data=sky_map)
