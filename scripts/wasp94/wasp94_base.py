# create an instrument with the appropriate settings
from mosasaurus.instruments import LDSS3C
i = LDSS3C(grism='vph-all')

# for working on a remote hardisk
import os
path = '/Users/zkbt/Cosmos/Data/Magellan/LDSS3/'

if os.path.exists(path):
    i.setupDirectories(path)

i.extractiondefaults['spatialsubarray'] = 200
i.extractiondefaults['narrowest'] = 4
i.extractiondefaults['widest'] = 20
i.extractiondefaults['numberofapertures'] = 5
i.extractiondefaults['zapcosmics'] = False
i.summarize()

# create a target, pulling values from simbad
from mosasaurus.Target import Target
import astropy.units as u
t = Target(starname='WASP-94A', name='WASP-94Ab')
t.summarize()
t.star.summarize()
