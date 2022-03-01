See the `README <https://github.com/christianhbye/lusee_sky_simulations/blob/main/README.rst>`_ for a general installation guide and how to get started. In general, all that is needed is to clone the repository and install the directory with pip.


Installing ARES
################

Models of the global 21-cm signal are generated with `ARES <https://github.com/mirochaj/ares>`_. This dependency is not automatically installed with the luseesky repository. Instead, the user should install ARES seperateley, following their `installation guide <https://github.com/mirochaj/ares#getting-started>`_. As of March 1, 2022, that involves cloning the repository and running the scripts setup.py and remote.py. For example:

.. code:: bash

    source .venv/bin/activate  # activate the virtual env used for luseesky
    git clone https://github.com/mirochaj/ares.git
    cd ares
    python setup.py install
    python remote.py


Installing ULSA
###############

If desired, `ULSA <https://github.com/Yanping-Cong/ULSA>`_ can be used to generate sky models described in `Cong et al, 2021 <https://ui.adsabs.harvard.edu/abs/2021ApJ...914..128C/abstract>`_. To make ULSA work, follow these steps:

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
