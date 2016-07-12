shortcuts = {}#{'/Users/zkbt/Cosmos/Data/Magellan/LDSS3':'...'}
import zachopy.Talker
zachopy.Talker.shortcuts = shortcuts
#zachopy.Talker.line = 200
Talker = zachopy.Talker.Talker

import astropy.io.fits, astropy.io.ascii, astropy.time
#from astropy.io import fits, ascii
#from astropy.time import Time

import matplotlib.pyplot as plt, numpy as np, matplotlib.animation
# ignore errors from divide by zero
np.seterr(divide='ignore')

#import zachopy.display, zachopy.oned
#from ds9 import *

import glob, os, string, copy, shutil

# ... to use
import scipy.interpolate, scipy.signal, scipy.integrate
#import scipy.interpolate as interp
#from scipy import interpolate
#from numpy.polynomial import polynomial as P

# ... to use utilities from the zachopy toolbox
import zachopy.utils, zachopy.regions, zachopy.iplot, zachopy.oned, zachopy.twod, zachopy.cmaps

# load tools to interface with ds9
#from zachopy.displays.ds9 import ds9 as zachods9
from zachopy.displays.loupe import loupe

from  zachopy.borrowed.mpfit.mpfit import mpfit

from Tools import *

mosasaurusdirectory = os.path.split(os.path.split(__file__)[0])[0] + '/'
