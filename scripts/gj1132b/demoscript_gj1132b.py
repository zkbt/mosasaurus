# create a target
from mosasaurus.Target import Target
t = Target(starname='GJ1132', name='GJ1132b')
t.summarize()

# create an instrument
from mosasaurus.instruments import LDSS3C
i = LDSS3C(grism='vph-red')
i.summarize()

# create a night
from mosasaurus.Night import Night
n = Night('ut160227_28', instrument=i)


# create an observation of this target, with this instrument, on this night
from mosasaurus.Observation import Observation
import numpy as np

o = Observation(t, i, n)
o.nReference = np.arange(1595, 1596+1)
o.nScience = np.arange(1625, 2025+1)
o.nHe = np.arange(162, 166+1)
o.nNe = np.arange(172, 176+1)
o.nAr = np.arange(177, 181+1)
o.nWideFlat = np.arange(194, 204+1)
o.nBias = np.arange(141, 151+1)
o.nWideMask = np.arange(156, 158+1)
o.nThinMask = np.arange(159, 161+1)
o.nDark = np.arange(205, 230+1)
o.nFinder = np.arange(1585, 1585+1)
