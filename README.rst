Tools for evaluating LuSEE-Night antenna models.

To install:

.. code:: bash

  python -m venv .venv  # make virtual env
  source .venv/bin/activate
  python -m pip install .  # use .[dev] to get extra developer dependencies

After running the above commands, do the following to make ULSA work if needed:

Download the Intel Fortran compiler from: https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran.
Then:

.. code:: bash

  cd ~/Downloads
  chmod u+rwx l_fortran-compiler_p_2022.0.2.83.sh  # makes the bash file exectuable
  ./l_fortran-compiler_p_2022.0.2.83.sh
 
This should open the Intel installer GUI. Follow the instructions.

Then, set the enviornment variables for :code:`ifort`:

.. code:: bash

  source ~/intel/oneapi/setvars.sh
  

Now, make the shared library (.so) file:

.. code:: bash

  cd .venv/lib/python3.8/site-packages/ULSA/NE2001
  unzip NE2001_4python.zip
  cd NE2001_4python/src.NE2001
  make so
  cd ../../../sky_map
  vim produce_absorbed_sky_map.py

Change line 39 to e.g.

.. code:: python

  libNE2001 = ct.CDLL('/home/christian/Documents/research/lusee/lusee_sky_simulations/.venv/lib/python3.8/site-packages/'
                      'ULSA/NE2001ULSA/NE2001_4python/src.NE2001/libNE2001.so')

Alternatively add the following lines:

.. code:: python

  from ULSA import __path__
  libNE2001 = ct.CDLL(__path__[0] + '/NE2001/NE2001/4_python/src.NE2001/libNE2001.so')

Finally, cd to lusee_sky_simulations and reinstall:

.. code:: bash

  python -m pip install .


To make a new skymap with ULSA, download the necessary data files from https://ulsa.readthedocs.io/en/latest/data_used.html and move them to ../obs_sky_data (relative to the ULSA installation path). Then copy all the files in ULSA/NE2001/NE2001_4python/bin_NE2001 to the directory you want to run ULSA in.
