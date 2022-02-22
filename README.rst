Tools for evaluating LuSEE-Night antenna models.

To install:

.. code:: bash

  python -m venv .venv  # make virtual env
  source .venv/bin/activate
  python -m pip install .  # use .[dev] to get extra developer dependencies

After running the above commands, do the following to make ULSA work:

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
  
Add the LD_LIBRARY_PATH environment variable to the Jupyter kernel for use in a Jupyter Notebook (https://stackoverflow.com/questions/37890898/how-to-set-env-variable-in-jupyter-notebook):

.. code:: bash

  jupyter kernelspec list  # find which kernel is related to the .venv virtual env. Might have to add one manually
  cd ~/.local/share/jupyter/kernels/.venv  # this is where my kernel is
  vim kernel.json
  
In the file kernel.json, add the key "env" with the value :code:`{"LD_LIBRARY_PATH": <library path>}` and set :code:`<library path>` to the desired value (I just copy-pasted the output of :code:`echo $LD_LIBRARY_PATH`).

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

