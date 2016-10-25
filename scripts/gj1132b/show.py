#!/usr/bin/env python

'''This script will compile 1D spectra into cubes, and visualize them.'''

# create (and save!) a cube, and visualize it
from mosasaurus.Cube import Cube
import sys

try:
    c = Cube(sys.argv[1])
    c.populate(remake=False, visualize=False)
    c.movieCube(stride=1)
except IndexError:
    print '''

    Example usage:
    ./show.py gj1132_0227.obs

    '''
