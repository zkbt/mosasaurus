from lhs1140_base import *

# create a night to analyze
from mosasaurus.Night import Night
n = Night('ut171026_27', instrument=i)
n.createNightlyLog(remake=False)

# create an observation
from mosasaurus.Observation import Observation
o = Observation(t, i, n)
o.setupFilePrefixes(science=['LHS1140 exposure'], reference=['LHS1140 Field'], flat=['lamp'])

# create a reducer to analyze this observation
from mosasaurus.Reducer import Reducer
r = Reducer(o, visualize=False)
r.reduce()

from mosasaurus.Cube import Cube
c = Cube(o, width=6)
c.setStars(target='aperture_546_1205', comparisons='aperture_165_1139')
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




#c.exportShiftStretch()

'''
c.imageCube(remake=True)
c.movieCube(stride=1, remake=True)
c.shiftCube()
c.imageCube(remake=True)
c.movieCube(stride=1, remake=True)
#c.nudgeWavelengths()
'''
