from datetime import datetime, timedelta
from pyuvsim import mpi, profiling, uvsim
import time

# Will add some methods to automatically write the yaml files

profiling.set_profiler(outfile_prefix="./profile", dump_raw=False)
t0 = time.time()
uvsim.run_uvsim("../sim_files/obsparam.yaml", return_uv=False, quiet=False)
dt = time.time() - t0
maxrss = mpi.get_max_mode_rss()
rtime = str(timedelta(seconds=dt))
if isinstance(maxrss, float):
    print("\tRuntime: {} \n\tMaxRSS: {:.3f} GiB".format(rtime, maxrss))
if hasattr(profiling.prof, "meta_file"):
    with open(profiling.prof.meta_file, "a") as afile:
        afile.write("Runtime \t {}\n MaxRSS \t {:.3f}\n".format(rtime, maxrss))
        afile.write("Date/Time \t {}".format(str(datetime.now())))
