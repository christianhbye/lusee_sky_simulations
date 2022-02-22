# import h5py
import numpy as np
from pyuvsim import uvsim

# Will add some methods to automatically write the yaml files

uvd = uvsim.run_uvsim("../sim_files/obsparam.yaml", return_uv=True, quiet=False)

max_imag = np.max(np.abs(uvd.data_array.imag))
print(f"{max_imag=}")

assert np.isclose(max_imag, 0.0, atol=1e-3)

uvd.data_array = uvd.data_array.real

save_dict = {
        "data": uvd.data_array,
        "freq": uvd.freq_array,
        "lst": uvd.lst_array
}

np.savez("./results", **save_dict)

#with h5py.File("./results.h5", "w") as hf:
#    hf.create_dataset("uvdata", data=save_dict)
