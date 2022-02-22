Tools for evaluating LuSEE-Night antenna models.

To install:

.. code:: bash

  python -m venv .venv  # make virtual env
  source .venv/bin/activate
  python -m pip install .  # use .[dev] to get extra developer dependencies

Note: One of the dependencies is ULSA which requires a specific installation. After running the above commands, do the following to make ULSA work:

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

