import sys
#sys.path.append('/home/hdiamond/local/lib/python2.7/site-packages/')
sys.path.append('/home/hdiamond/anaconda3/lib/python3.6/site-packages/')
sys.path.append('/h/mulan0/code/')
sys.path.append('/h/mulan0/code/mosasaurus')
#sys.path.append('/h/mulan0/code/mosasaurus/mosasaurus/craftroom')


# create an instrument with the appropriate settings
from mosasaurus.instruments import LDSS3C
i = LDSS3C(grism='VPH-red')
#from mosasaurus.instruments import IMACS
#i = IMACS(grism='Gri-300-26.7')

# for working on a remote hardisk
import os
path = '/h/mulan0/data/LHS3844/LDSS3C/' # directory that contains /data/

if os.path.exists(path):
    i.setupDirectories(path)

i.extractiondefaults['spatialsubarray'] = 75
i.extractiondefaults['narrowest'] = 6
i.extractiondefaults['widest'] = 10
i.extractiondefaults['numberofapertures'] = 3
i.extractiondefaults['zapcosmics'] = False
i.summarize()

# create a target, pulling values from simbad
from mosasaurus.Target import Target
import astropy.units as u
t = Target(starname='LHS3844', name='LHS3844b')
t.summarize()
t.star.summarize()