import ares  # type: ignore
import numpy as np

# params to save:
blobs_scalar = ['z_D']  # redshift
blobs_1d = ['dTb']  # brightness temp
freqs = np.linspace(1e6, 50e6, num=50)
blobs_1d_z = 1420e6/freqs-1  # redshifts

base_pars = {
 'problem_type': 101,
 'tanh_model': False,
 'blob_names': [blobs_scalar, blobs_1d],
 'blob_ivars': [None, [('z', blobs_1d_z)]],
 'blob_funcs': None,
}

# refvals = {'Nlw': 1e4, 'Nion': 4e3, 'fX': 0.2, 'fstar': 0.01}

mg = ares.inference.ModelGrid(**base_pars)

axes = {
      #  'Nlw': np.array([1e1, 1e3, 1e4, 1e5, 1e6]),
      #  'Nion': np.hstack([np.array([1e3, .5e4]), np.linspace(1e4, 1e5, 6)]),
        'fX': np.logspace(-2, 3, 8), 
        'fstar': np.logspace(np.log10(5e-3), np.log10(5e-1), 5),
        'Tmin': np.logspace(4, np.log10(5e6), 8)
    }

mg.axes = axes

if __name__=="__main__":
    mg.run("test_2d_grid", clobber=True)
