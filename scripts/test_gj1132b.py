
# create an instrument with the appropriate settings
from mosasaurus.instruments import LDSS3C
i = LDSS3C(grism='vph-red')
i.extractiondefaults['zapcosmics'] = False
i.summarize()

# create a night to analyze
from mosasaurus.Night import Night
n = Night('ut160227_28', instrument=i)
n.createNightlyLog(remake=False)

# create a target from values (in case you want to work offline)
from mosasaurus.Target import Target
import astropy.units as u
t = Target(starname='GJ1132', name='GJ1132b',
            ra=153.71622708*u.deg, dec=-47.15684639*u.deg,
            pmra=-1046., pmdec=416., epoch=2000.0)
t.summarize()
t.star.summarize()

# create an observation
from mosasaurus.Observation import Observation
o = Observation(t, i, n)

# create a reducer to analyze this observation
from mosasaurus.Reducer import Reducer
r = Reducer(o, visualize=False)
r.visualize = True
r.reduce()

from mosasaurus.Cube import Cube
c = Cube(o, width=6)
c.populate(shift=False)
c.setStars(target='aperture_100_809', comparisons=['aperture_693_656'])
c.shiftCube(plot=True)
#c.nudgeWavelengths()
#c.movieCube(stride=1, figsize=(12,8))
#c.listWidths()
#c.populate(shift=False)
#
