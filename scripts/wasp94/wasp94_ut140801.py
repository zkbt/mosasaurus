from wasp94_base import *

# create a night to analyze
from mosasaurus.Night import Night
n = Night('ut140801', instrument=i)
n.createNightlyLog(remake=False)

# create an observation
from mosasaurus.Observation import Observation
o = Observation(t, i, n)
o.setupFilePrefixes(science=['WASP-94'], reference=['WASP-94'], flat=['flat'])

# create a reducer to analyze this observation
from mosasaurus.Reducer import Reducer
r = Reducer(o, visualize=False)
r.reduce()

from mosasaurus.Cube import Cube
c = Cube(o, width=16)
c.populate(shift=False, max=None)
c.setStars(target='aperture_709_1066', comparisons='aperture_751_1066')
c.save()

from mosasaurus.WavelengthRecalibrator import WavelengthRecalibrator
wr = WavelengthRecalibrator(c)

#c.imageCube(keys=['raw_counts'], stars=[c.target])
#c.imageCube()
#c.populate(shift=True, max=None)
#c.imageCube(keys=['raw_counts'], stars=[c.target])
#c.imageCube()
#c.exportShiftStretch()



#c.shiftCube()
#c.imageCube(keys=['raw_counts'], stars=[c.target], remake=True)


'''
c.movieCube(stride=1, remake=False)
c.imageCube(remake=True)
c.movieCube(stride=1, remake=False)
'''
#c.nudgeWavelengths()
