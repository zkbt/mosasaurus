#!/usr/bin/env python


'''This script will reduce LDSS3 raw data into 1D spectra.'''


from mosasaurus.Reducer import Reducer
import sys
try:
    r = Reducer(sys.argv[1])
    r.reduce()
except IndexError:
    print '''

Example usage:
    ./reduce.py gj1132_0227.obs

'''
