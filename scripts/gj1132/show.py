#!/usr/bin/env python

'''This script will compile 1D spectra into cubes, and visualize them.'''

# create (and save!) a cube, and visualize it
import sys
sys.path.append('/home/hdiamond/local/lib/python2.7/site-packages/')
sys.path.append('/h/mulan0/code/')
sys.path.append('/h/mulan0/code/mosasaurus')
from mosasaurus.Cube import Cube


try:
    c = Cube(sys.argv[1])
    c.populate(remake=False, visualize=False, shift=False)
    c.movieCube(stride=1)
except IndexError:
    print('''

    Example usage:
    ./show.py gj1132_0227.obs

    ''')
