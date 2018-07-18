from lhs1140_base import *

# create a night to analyze
from mosasaurus.Night import Night
n = Night('ut171026_27', instrument=i)
n.createNightlyLog(remake=False)

# create an observation
from mosasaurus.Observation import Observation
o = Observation(t, i, n)
o.setupFilePrefixes(science=['LHS1140'], reference=['LHS1140'], flat=['lamp'])

# create a reducer to analyze this observation
from mosasaurus.Reducer import Reducer
r = Reducer(o, visualize=False)
r.reduce()

from mosasaurus.Cube import Cube
c = Cube(o, width=6)
c.setStars(target='aperture_1455_3430', comparisons=['aperture_795_1880', 'aperture_1352_1515', 'aperture_1811_3491', 'aperture_3808_2678'])
c.populate(shift=False, max=None)
c.savable=c.savable + ['target', 'comparisons']
c.save()
c.imageCube(keys=['raw_counts'], stars=[c.target])
#c.imageCube()
#c.populate(shift=True, max=None)
#c.imageCube(keys=['raw_counts'], stars=[c.target])
#c.save()
#c.imageCube()
#c.exportShiftStretch()

from mosasaurus.WavelengthRecalibrator import WavelengthRecalibrator
wr = WavelengthRecalibrator(c, visualize=True)

# fix up the wavelength calibration for each exposure
r.mask.setup()
r.mask.addWavelengthCalibration(shift=True)

# repopulate the cube
c.populate(shift=True, remake=True)
c.save()
c.imageCube(keys=['raw_counts'], stars=[c.target])

#c.exportShiftStretch()

'''
c.imageCube(remake=True)
c.movieCube(stride=1, remake=True)
c.shiftCube()
c.imageCube(remake=True)
c.movieCube(stride=1, remake=True)
#c.nudgeWavelengths()
'''
