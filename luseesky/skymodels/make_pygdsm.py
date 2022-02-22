from astropy.coordinates import Galactic
from astropy_healpix import HEALPix
from astropy import units as u
import numpy as np
from pygdsm import GlobalSkyModel2016
import pyradiosky

print("Generate GSM.")
REF_FREQ = 25  # MHz
gsm = GlobalSkyModel2016(freq_unit='MHz')
gsm.generate(REF_FREQ)

NSIDE = int(np.sqrt(gsm.generated_map_data.size/12))
NPIX = gsm.generated_map_data.size
HPX_INDS = np.arange(NPIX)
hp = HEALPix(nside=NSIDE, order="ring", frame=Galactic())
assert hp.npix == NPIX
coords = hp.healpix_to_skycoord(HPX_INDS).transform_to("icrs")
ra, dec = coords.ra, coords.dec
# some mock stokes params assuming Ex == Ey
stokes = u.Quantity(np.zeros((4, 1, NPIX)), unit=u.K)
stokes[0, 0] = gsm.generated_map_data * u.K  # set I to the GDSM Temp
stokes[2, 0] = stokes[0, 0]  # U = I, Q=V=0
print("Convert to pyradiosky SkyModel.")
skymodel = pyradiosky.SkyModel()
skymodel.Ncomponents = NPIX
skymodel.Nfreqs = 1
skymodel.coherency_radec = pyradiosky.utils.stokes_to_coherency(stokes)
skymodel.component_type = "healpix"
skymodel.dec = dec
skymodel.history = "PyGDSM2016"
skymodel.ra = ra
skymodel.spectral_type = "spectral_index"
skymodel.stokes = stokes

skymodel.hpx_inds = HPX_INDS
skymodel.hpx_order = "ring"
skymodel.nside = NSIDE
skymodel.reference_frequency = REF_FREQ * np.ones(NPIX) * u.MHz
skymodel.spectral_index = -2.5 * np.ones(NPIX)
assert skymodel.check()

print("Convert to effective point sources.")
skymodel.healpix_to_point()
print("Filter to just the brightest sources.")
# the brightest sources has greatest stokes I
comp_ind = np.argsort(stokes[0, 0])[-50:]  # 100 brightest
skymodel.select(comp_ind, inplace=True)
print("Save catalog to text.")
skymodel.write_text_catalog("./pygdsm16_srcs.txt")
