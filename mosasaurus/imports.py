shortcuts = {}#{'/Users/zkbt/Cosmos/Data/Magellan/LDSS3':'...'}
import craftroom.Talker
craftroom.Talker.shortcuts = shortcuts
#craftroom.Talker.line = 200
Talker = craftroom.Talker.Talker

import astropy.io.fits, astropy.io.ascii, astropy.time
from astropy import time, coordinates as coord, units as u

import matplotlib.pyplot as plt, numpy as np, matplotlib.animation
plt.switch_backend('qt5Agg')

# ignore errors from divide by zero
np.seterr(divide='ignore')

#import craftroom.display, craftroom.oned
#from ds9 import *

import glob, os, string, copy, shutil, warnings

# ... to use
import scipy.interpolate, scipy.signal, scipy.integrate
#import scipy.interpolate as interp
#from scipy import interpolate
#from numpy.polynomial import polynomial as P

# ... to use utilities from the craftroom toolbox
import craftroom.displays.regions, craftroom.displays.iplot, craftroom.oned, craftroom.twod, craftroom.cmaps, craftroom.resample
from craftroom.resample import fluxconservingresample

# load tools to interface with ds9
#from craftroom.displays.ds9 import ds9 as zachods9
from craftroom.displays.loupe import loupe

#from craftroom.borrowed.mpfit.mpfit import mpfit

from tqdm import tqdm

from .Tools import *

mosasaurusdirectory = os.path.split(os.path.split(__file__)[0])[0] + '/'
