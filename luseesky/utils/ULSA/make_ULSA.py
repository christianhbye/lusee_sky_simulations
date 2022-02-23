# type: ignore

from astropy.io import fits
import numpy as np
from ULSA.sky_map.produce_absorbed_sky_map import absorption_JRZ

FREQ = 25
NSIDE = 32
INDEX_TYPE = "constant_index"

# default params
DISTANCE = 50  # kpc
PAR_408 = np.array(
    [43.09932087, 3.40820378, 0.46230938, 1.11894379, 1.2344227]
)
RAW_DIFFUSE = False
FDIRS = None
DEF_PARS = True
IN_SPEC_IDX = None
CRIT_DI = False
OUT_FREE = False

sky_map = absorption_JRZ(
    v=FREQ,
    nside=NSIDE,
    index_type=INDEX_TYPE,
    distance=DISTANCE,
    using_raw_diffuse=RAW_DIFFUSE,
    v_file_dir=FDIRS,
    using_default_params=DEF_PARS,
    input_spectral_index=IN_SPEC_IDX,
    params_408=PAR_408,
    critical_dis=CRIT_DI,
    output_absorp_free_skymap=OUT_FREE,
)

sky_map.mpi()

save_params = {
    "freqs_MHz": FREQ,
    "nside": NSIDE,
    "index_type": INDEX_TYPE,
    "distance_kpc": DISTANCE,
    "params_408": PAR_408,
    "raw_diffuse": RAW_DIFFUSE,
    "freq_file_dirs": FDIRS,
    "default_params": DEF_PARS,
    "input_spectral_idx": IN_SPEC_IDX,
    "crtiical_dist": CRIT_DI,
    "output_abs_free": OUT_FREE,
}


hdu = fits.PrimaryHDU(sky_map)
hdu.writeto("../../skymodel.skymap.fits", header=save_params)
