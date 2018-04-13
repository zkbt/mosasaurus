import sys
#sys.path.append('/home/hdiamond/local/lib/python2.7/site-packages/')
sys.path.append('/home/hdiamond/anaconda3/lib/python3.6/site-packages/')
sys.path.append('/h/mulan0/code/')
sys.path.append('/h/mulan0/code/mosasaurus')
#sys.path.append('/h/mulan0/code/mosasaurus/mosasaurus/craftroom')


# create an instrument with the appropriate settings
from mosasaurus.instruments import LDSS3C
i = LDSS3C(grism='VPH-red')

# for working on a remote hardisk
import os
path = '/h/mulan0/data/LHS1140/LDSS3C/' # directory that contains /data/

if os.path.exists(path):
    i.setupDirectories(path)

i.extractiondefaults['spatialsubarray'] = 50
i.extractiondefaults['narrowest'] = 4
i.extractiondefaults['widest'] = 12
i.extractiondefaults['numberofapertures'] = 5
i.extractiondefaults['zapcosmics'] = False
#del i.alignmentranges['Ca1'] # to remove this key
i.summarize()

# create a target, pulling values from simbad
from mosasaurus.Target import Target
import astropy.units as u
t = Target(starname='LHS1140', name='LHS1140b')
t.summarize()
t.star.summarize()
