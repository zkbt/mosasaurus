'''
Tools to organize data downloaded from MAST

'''
import os, glob
import numpy as np
from craftroom.utils import mkdir

def create_visit_directories(path='HST/'):
    '''
    Organize files from the directories downloaded
    using a curl script through the MAST portal
    (with files inside a directory for each exposure)
    into unique visits.

    (This are different from the staged downloads,
    in which all files are in a single directory.)
    '''

    # find all the directories
    directories = glob.glob(os.path.join(path, '*'))

    # pull out the 6-letter visit prefix
    prefixes = [d.split('/')[1][0:6] for d in directories]

    # identify the unique visits
    unique = np.unique(prefixes)
    visits = {f'visit{i+1}':p for i, p in enumerate(unique)}

    # create a directory for each visit and populate it
    for v in visits:
        mkdir(v)
        thisvisit = glob.glob(f'HST/{visits[v]}*')
        for t in thisvisit:
            print(t, '->', f'{v}/{e}')
            os.rename(t, f'{v}/{e}')
            e = t.split('/')[-1]
