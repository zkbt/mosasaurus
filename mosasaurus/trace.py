

import numpy as np
import mosasaurus.gaussian as g
import scipy.interpolate as spi
import matplotlib.pyplot as plt
import mosasaurus.weighted as weighted

def calc_trace(data, mask, nknots=50, deg=3, guess=None, isplots=0):
    '''
    Calculate the trace of a spectrum.  Assumes y is dispersion direction and x is spatial direction.
    
    PARAMETERS
    ----------
    data    : 2D array
    mask    : 2D array
    nknots  : Number of knots in Gaussian fit
    deg     : Degree polynomial for spline interpolation
    guess   : starting guess position for trace
    
    RETURNS
    -------
    trace   : 1D array
    
    '''
    ny,nx = data.shape
    centers = np.zeros(nknots)
    knots   = np.zeros(nknots)
    for i in range(nknots):
        istart  = i*ny/nknots
        iend    = (i+1)*ny/nknots
        #subdata = np.sum(data[istart:iend+1]*mask[istart:iend+1],axis=0)
        subdata = weighted.mean(data[istart:iend+1],1./mask[istart:iend+1],axis=0)
        #subdata = np.sum(data[istart:iend+1]*(1./mask[istart:iend+1])**2, 0)/np.sum((1./mask[istart:iend+1])**2, 0)
        subdata[np.where(np.isnan(subdata))] = 0
        submask = np.ones(subdata.shape, dtype=int)
        submask[np.where(subdata == 0)] = 0
        if guess == None:
            guess   = np.argmax(subdata)
        if submask.sum() > 0:
            p, err  = g.fitgaussian(subdata,np.arange(nx),guess=(nx/10.,int(guess),subdata.max()),mask=submask)
            centers[i]  = guess = p[1]
        else:
            centers[i]  = guess
        knots[i]    = (istart+iend)/2
        if isplots == 4:
            plt.figure(1)
            plt.clf()
            plt.suptitle(str(i))
            plt.plot(np.where(submask)[0],subdata[np.where(submask)],'bo')
            ymin,ymax=plt.ylim()
            plt.vlines(centers[i], ymin, ymax)
            plt.pause(0.1)
        #plt.plot(subdata)
    
    spl     = spi.UnivariateSpline(knots, centers, k=deg, s=0)
    trace   = spl(range(ny))
    
    return trace

def realign(subdata, submask, variance, bg, spec_width, trace):
    '''

    Parameters
    ----------
    subdata     : background-subtracted data subarray (of dimensions subny, subnx)
    submask     : mask subarray, 1 = good, 0 = bad
    variance    : uncertainty^2 subarray
    bg          : background subarray
    trace       : trace array of size subny
    spec_width  : half-width of aperture size to use in optimal spectral extraction

    Returnsb
    -------
    trdata      : interpolated (trace-corrected) data array (of dimensions subny, spec_width*2)
                  used as input to optimal spectral extraction routine
    trmask      : interpolated mask array
    trvar       : interpolated variance array
    trbg        : interpolated background array
    '''

    subny,subnx = subdata.shape

    trdata  = np.zeros((subny,spec_width*2))
    trmask  = np.zeros((subny,spec_width*2),dtype=bool)
    trvar   = np.zeros((subny,spec_width*2))
    trbg    = np.zeros((subny,spec_width*2))
    spldata = spi.RectBivariateSpline(range(subny),range(subnx),subdata, kx=3, ky=3, s=0)
    splmask = spi.RectBivariateSpline(range(subny),range(subnx),submask, kx=3, ky=3, s=0)
    splvar  = spi.RectBivariateSpline(range(subny),range(subnx),variance, kx=3, ky=3, s=0)
    splbg   = spi.RectBivariateSpline(range(subny),range(subnx),bg, kx=3, ky=3, s=0)
    for i in range(subny):
        xi          = (np.arange(trace[i]-spec_width,trace[i]+spec_width)+0.5)[:spec_width*2]
        yi          = np.zeros(len(xi))+i
        trdata[i]   = spldata.ev(yi,xi)
        trmask[i]   = np.round(splmask.ev(yi,xi),0).astype(bool)
        trvar [i]   = splvar .ev(yi,xi)
        trbg  [i]   = splbg  .ev(yi,xi) 
    
    return trdata, trmask, trvar, trbg

