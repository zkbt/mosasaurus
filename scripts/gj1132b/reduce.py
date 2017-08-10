#!/usr/bin/env pythonw


'''This script will reduce LDSS3 raw data into 1D spectra.'''


import sys
sys.path.append('/home/hdiamond/local/lib/python2.7/site-packages/')
sys.path.append('/h/mulan0/code/')
sys.path.append('/h/mulan0/code/mosasaurus')
from mosasaurus.Reducer import Reducer

try:
    r = Reducer(sys.argv[1])
    r.reduce()
except IndexError:
    print '''

Example usage:
    ./reduce.py gj1132_0227.obs

'''
