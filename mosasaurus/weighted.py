
import numpy as np

def mean(val, err, axis=None):
    '''
    Compute weighted mean.
    
    Parameters
    ----------
    val : Array of values
    err : Array of uncertainties
    '''
    val = np.array(val)
    err = np.array(err)
    return np.sum(val/err**2,axis=axis)/np.sum(1/err**2,axis=axis)

def error(err, axis=None):
    '''
    Compute standard error.
    '''
    err = np.array(err)
    return np.sqrt(1./np.sum(1./err**2,axis=axis))
