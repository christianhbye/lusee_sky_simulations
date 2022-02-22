Tools for evaluating LuSEE-Night antenna models.

To install:
python -m venv .venv  # make virtual env
source .venv/bin/activate
python -m pip install .  # use .[dev] to get extra developer dependencies

Note: One of the dependencies is ULSA which requires a specific installation. After running the above commands, do the following to make ULSA work:

cd .venv/lib/python3.8/site-packages/ULSA/NE2001
unzip NE2001\_4python.zip
cd NE2001\_4python/src.NE2001
make so
cd ../../../sky\_map
vim produce\_absorbed\_sky\_map.py

Change line 39 to e.g.
libNE2001 = ct.CDLL('/home/christian/Documents/research/lusee/lusee\_sky\_simulations/.venv/lib/python3.8/site-packages/ULSA/NE2001ULSA/NE2001\_4python/src.NE2001/libNE2001.so')

Alternatively add the following lines:
from ULSA import __path__
libNE2001 = ct.CDLL(__path__[0] + '/NE2001/NE2001/4\_python/src.NE2001/libNE2001.so')

Finally, cd to lusee\_sky\_simulations and reinstall:
python -m pip install .

