from dataclasses import dataclass, field
import fitsio
import numpy.typing as npt

@dataclass
class Beam:
    fname: str
    frequencies: npt.ArrayLike = field(init=False)
    fits = fitsio.FITS(se
