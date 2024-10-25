shortcuts = {}  # {'/Users/zkbt/Cosmos/Data/Magellan/LDSS3':'...'}
from .craftroom.Talker import Talker

Talker.shortcuts = shortcuts
# craftroom.Talker.line = 200
# Talker = craftroom.Talker.Talker

import astropy.io.fits, astropy.io.ascii, astropy.time
from astropy import time, coordinates as coord, units as u

import matplotlib.pyplot as plt, numpy as np, matplotlib.animation

# plt.switch_backend('qt5Agg')

# ignore errors from divide by zero
np.seterr(divide="ignore")

# import craftroom.display, oned
# from ds9 import *

import glob, os, string, copy, shutil, warnings

# ... to use
import scipy.interpolate, scipy.signal, scipy.integrate

# import scipy.interpolate as interp
# from scipy import interpolate
# from numpy.polynomial import polynomial as P

# ... to use utilities from the craftroom toolbox
from .craftroom.displays.iplot import iplot
from .craftroom import oned, twod


from chromatic import one2another
from chromatic.resampling import resample_while_conserving_flux


# load tools to interface with ds9
# from .craftroom.displays.ds9 import ds9 as zachods9
from .craftroom.displays.loupe import loupe

# from .craftroom.borrowed.mpfit.mpfit import mpfit

from tqdm import tqdm

from .Tools import *

mosasaurusdirectory = os.path.split(os.path.split(__file__)[0])[0] + "/"


def clean(s):
    cleaned = s + ""
    nasty = " !@#$%^&*()"
    for character in nasty:
        cleaned = cleaned.replace(character, "")
    return s
