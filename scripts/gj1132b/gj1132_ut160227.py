from gj1132_base import *

# create a night to analyze
from mosasaurus.Night import Night
n = Night('ut160227_28', instrument=i)
n.createNightlyLog(remake=False)

# create an observation
from mosasaurus.Observation import Observation
o = Observation(t, i, n)
o.setupFilePrefixes(science=['GJ1132 science'], reference=['GJ1132 Field'], flat=['lamp'])

# create a reducer to analyze this observation
from mosasaurus.Reducer import Reducer
r = Reducer(o, visualize=False)
r.reduce()

from mosasaurus.Cube import Cube
c = Cube(o, width=16)
c.setStars(target='aperture_99_809', comparisons=['aperture_181_1018', 'aperture_542_1024', 'aperture_692_654', 'aperture_832_1173'])
c.populate(shift=False, max=None)
c.savable=c.savable + ['target', 'comparisons']
c.save()
#c.imageCube(keys=['raw_counts'], stars=[c.target])
#c.imageCube()
#c.populate(shift=True, max=None)
#c.imageCube(keys=['raw_counts'], stars=[c.target])
#c.save()
#c.imageCube()
#c.exportShiftStretch()

from mosasaurus.WavelengthRecalibrator import WavelengthRecalibrator
wr = WavelengthRecalibrator(c)

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
