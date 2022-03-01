from astropy.coordinates import Galactic
from astropy_healpix import HEALPix
from astropy import units as u
import healpy
import numpy as np
from pygdsm import GlobalSkyModel2016
import pyradiosky
from typing import Optional
import warnings


def nside2npix(nside: int) -> int:
    return 12*nside**2

def npix2nside(npix: int) -> int:
    nside = np.sqrt(npix/12)
    notint = ~np.isclose(nside, int(nside))
    notpower = ~np.isclose(np.log2(nside), int(np.log2(nside)))
    if notint or notpower:
        raise ValueError(f"npix has an invalid value {npix}.")
    return int(nside)

def _degrade(sky_map: np.ndarray, nside_out: int) -> np.ndarray:
    return healpy.ud_grade(
            sky_map,
            nside_out=nside_out,
            order_in="RING",
            order_out="RING"
            )

class Skymodel:
    
    def __init__(
            frequency: float = 25,
            freq_unit: str = "MHz",
            nside: Optional[int] = 16,
            npix: Optional[int] = None,
            healpix_map: Optional[np.ndarray] = None,
            degrade: bool = False,
            base_fname: str = "./pygdsm16"
            ):
        """
        Must at least specify one of nside, npix, or healpix_map. If they
        are incompatible, nside and npix will change.

        nside: int, nside of healpix map, must be power of 2
        npix: int, npix of healpix map, must be s.t. sqrt(npix/12)
        is power of 2
        healpix_map: np.ndarray, a sky map in RING healpix coords
        degrade: bool, whether to degrade the healpix map to the given nside
        """
     
        if self.healpix_map is not None:
            if self.nside is None and self.npix is not None:
                self.nside = npix2nside(self.npix)
            if self.degrade:
                self.healpix_map = _degrade(self.healpix_map, self.nside)
            self.degrade = False
            self.npix = healpix_map.size
            self.nside = npix2nside(npix)
        elif self.nside is None:
            if self.npix is None:
                raise ValueError(
                        "Must specify at least one of nside, npix, and "
                         "healpix_map."
                         )
            else:
                self.nside = npix2nside(self.npix)
        elif self.npix is None:
            self.npix = nside2npix(self.nside)

        freq_conversion = {"Hz": 1e-6, "kHz": 1e-3, "MHz:" 1., "GHz:" 1e3}
        if not self.freq_unit in freq_conversion:
            raise ValueError(
            f"Invalid frequency unit {freq_unit}, must be in"
            f"{list[freq_conversion.keys()]}"
            )
        self.frequency *= freq_conversion[self.freq_unit]
        self.freq_unit = "MHz"

    def gen_gsm(self):
        gsm = GlobalSkyModel2016(freq_unit="MHz")
        gsm.generate(frequency)
        healpix_map = gsm.generated_map_data
        # degrading map
        if degrade:
            healpix_map = _degrade(healpix_map, self.nside)
        self.healpix_map = healpix_map

    def make_pyradio_skymap(self):
        HPX_INDS = np.arange(self.npix)
        hp = HEALPix(nside=self.nside, order="ring", frame=Galactic())
        coords = hp.healpix_to_skycoord(HPX_INDS).transform_to("icrs")
        ra, dec = coords.ra, coords.dec
        # some mock stokes params assuming Ex == Ey
        stokes = u.Quantity(np.zeros((4, 1, NPIX)), unit=u.K)
        stokes[0, 0] = dgr_map * u.K  # set I to the GDSM Temp
        # stokes[2, 0] = stokes[0, 0]  # U = I, Q=V=0
        skymodel = pyradiosky.SkyModel()
        skymodel.Ncomponents = self.npix
        skymodel.Nfreqs = 1
        skymodel.coherency_radec = pyradiosky.utils.stokes_to_coherency(stokes)
        skymodel.component_type = "healpix"
        skymodel.dec = dec
        skymodel.history = "Degraded PyGDSM2016"
        skymodel.ra = ra
        skymodel.spectral_type = "spectral_index"  # XXX: make flexible
        skymodel.stokes = stokes
        skymodel.hpx_inds = HPX_INDS
        skymodel.hpx_order = "ring"
        skymodel.nside = self.nside
        skymodel.reference_frequency = self.frequency * np.ones(self.npix) * u.MHz
        skymodel.spectral_index = -2.5 * np.ones(self.npix)  # XXX
        assert skymodel.check()
        skymodel.healpix_to_point()  # convert to point sources
        skymodel.write_text_catalog(self.base_name + f"_nside{NSIDE}.txt")

if __name__=="__main__":
    sm = SkyModel()
    if sm.healpix_map is None:
        sm.gen_gsm()
    sm.make_pyradio_skymap()
