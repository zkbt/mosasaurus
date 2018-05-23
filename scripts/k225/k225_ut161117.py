from k225_base import *

# create a night to analyze
from mosasaurus.Night import Night
n = Night('ut161117_18', instrument=i)
n.createNightlyLog(remake=False)

# create an observation
from mosasaurus.Observation import Observation
o = Observation(t, i, n)
o.setupFilePrefixes(science=['K225'], reference=['K225'], flat=['lamp'])

# create a reducer to analyze this observation
from mosasaurus.Reducer import Reducer
r = Reducer(o, visualize=False)
r.reduce()

from mosasaurus.Cube import Cube
c = Cube(o, width=6)
c.setStars(target='aperture_1629_2762', comparisons=['aperture_218_1134', 'aperture_610_2849', 'aperture_1216_2293', 'aperture_1716_1553', 'aperture_1999_1263', 'aperture_2625_2897' 'aperture_2922_2427', 'aperture_2937_797', 'aperture_3003_2953', 'aperture_3256_2992', 'aperture_3602_2199', 'aperture_3801_3277'])
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
wr = WavelengthRecalibrator(c)

# fix up the wavelength calibration for each exposure
r.mask.setup()
r.mask.addWavelengthCalibration(shift=True)

# repopulate the cube
c.populate(shift=True, remake=True)
c.save()
c.imageCube(keys=['raw_counts'], stars=[c.target])




'''
c.imageCube(remake=True)
c.movieCube(stride=1, remake=True)
c.shiftCube()
c.imageCube(remake=True)
c.movieCube(stride=1, remake=True)
#c.nudgeWavelengths()
'''
