__version__ = '0.2'

# specify whether we're calling this from within setup.py
try:
    __MOSASAURUS_SETUP__
except NameError:
    __MOSASAURUS_SETUP__ = False

if not __MOSASAURUS_SETUP__:
    # (run this stuff if it's not form within setup.py)
    pass
