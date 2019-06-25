
# create an instrument with the appropriate settings
from mosasaurus.instruments import DIS
i = DIS(grating='B400')
i.extractiondefaults['zapcosmics'] = False
i.summarize()

# create a night to analyze
from mosasaurus.Night import Night
n = Night('UT161114', instrument=i)
n.createNightlyLog(remake=False)

# create a target from values (in case you want to work offline)
from mosasaurus.Target import Target
import astropy.units as u
t = Target(starname='Kepler-51', name='Kepler-51')

t.summarize()
t.star.summarize()

# create an observation
from mosasaurus.Observation import Observation
o = Observation(t, i, n)
# FIXME! still need to make the directories encode the R/B cameras

# create a reducer to analyze this observation
from mosasaurus.Reducer import Reducer
r = Reducer(o, visualize=False)
r.visualize = True
r.reduce()

'''
from mosasaurus.Cube import Cube
c = Cube(o, width=6)
c.populate(shift=False)
c.movieCube(stride=1, remake=False)
c.setStars(target='aperture_78_0', comparisons=['aperture_390_0'])
c.save()
c.imageCube(keys=['raw_counts'], stars=[c.target])



from mosasaurus.WavelengthRecalibrator import WavelengthRecalibrator
wr = WavelengthRecalibrator(c, visualize=True)

# fix up the wavelength calibration for each exposure
r.mask.setup()
r.mask.addWavelengthCalibration(shift=True)

# repopulate the cube
c.populate(shift=True, remake=True)
c.imageCube(keys=['raw_counts'], stars=[c.target])
c.save()

# make movie of the shifted cube
c.movieCube(stride=1, remake=False)
'''
