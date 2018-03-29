# $Author: patricio $
# $Revision: 285 $
# $Date: 2010-06-18 17:59:25 -0400 (Fri, 18 Jun 2010) $
# $HeadURL: file:///home/esp01/svn/code/python/branches/patricio/photpipe/lib/gaussian.py $
# $Id: gaussian.py 285 2010-06-18 21:59:25Z patricio $

#! /usr/bin/env python

'''
    Name
    ----
    gaussian

    File
    ----
    gaussian.py

    Description
    -----------
    Routines for evaluating, estimating parameters of, and fitting Gaussians.

    Package Contents
    ----------------
    N-dimensional functions:

    gaussian(x, width=1., center=0., height=None, params=None)
        Evaluate the Gaussian function with given parameters at x
        (n-dimensional).
    fitgaussian(y, x)
        Calculates a Gaussian fit to (y, x) data, returns (width,
        center, height).

    1-dimensional functions:
    
    gaussianguess(y, x=None)
        Crudely estimates the parameters of a Gaussian that fits the
        (y, x) data.

    Examples:
    ---------
    See fitgaussian() example.

    Revisions
    ---------
    2007-09-17 0.1 jh@physics.ucf.edu Initial version 0.01, portions
                   adapted from http://www.scipy.org/Cookbook/FittingData.
    2007-10-02 0.2 jh@physics.ucf.edu Started making N-dimensional,
                   put width before center in args.
    2007-11-13 0.3 jh@physics.ucf.edu Made N-dimensional.
    2008-12-02 0.4 nlust@physics.ucf.edu Made fit gaussian return errors, and
                   fixed a bug generating initial guesses
    2009-10-25 0.5 jh@physics.ucf.edu Standardized all headers, fixed
                   an error in a fitgaussian example, added example
                   ">>>"s and plot labels.
'''

import numpy as np
import scipy.optimize as so
import disk as d

def gaussian(x, width=1.0, center=0.0, height=None, bgpars=[0.0, 0.0, 0.0]):
  """
    Evaluates the Gaussian and a background with given parameters at
    locations in x.

    Parameters
    ----------
    x : ndarray (any shape)
        Abcissa values.  Arranged as the output of np.indices() but
        may be float.  The highest dimension must be equal to the
        number of other dimensions (i.e., if x has 6 dimensions, the
        highest dimension must have length 5, and each of those must
        give the coordinate along the respective axis).  May also be
        1-dimensional.  Default: np.indices(y.shape).

    width : array_like
        The width of the Gaussian function, sometimes called sigma.
        If scalar, assumed constant for all dimensions.  If array,
        must be linear and the same length as the first dimension of
        x.  In this case, each element gives the width of the function
        in the corresponding dimension.  Default: [1.].

    center : array_like
        The mean value of the Gaussian function, sometimes called x0.
        Same scalar/array behavior as width.  Default: [0.].

    height : scalar
        The height of the Gaussian at its center.  If not set,
        initialized to the value that makes the Gaussian integrate to
        1.  If you want it to integrate to another number, leave
        height alone and multiply the result by that other number
        instead.  Must be scalar.  Default: [product(1./sqrt(2 * pi *
        width**2))].

    bgpars : ndarray or tuple, 3-element
        Background parameters, the elements determine a X- and Y-linearly
        dependant level, of the form:
        f = Y*bgparam[0] + X*bgparam[1] + bgparam[2]
        (Not tested for 1D yet).

    Returns
    -------
    results : ndarray, same shape as x (or first element of x if
        multidimensional)
        This function returns the Gaussian function of the given
        width(s), center(s), and height applied to its input plus a
        linear background level.  The Gaussian function is: f(x) =
        1./sqrt(2 * pi * width**2) * exp(-0.5 * ((x - center) /
        width)**2).  It is defined in multiple dimensions as the
        product of orthogonal, single-dimension Gaussians.

    Examples
    --------
    
    >>> import matplotlib.pyplot as plt
    >>> import gaussian as g
    
    >>> x = np.arange(-10., 10.005, 0.01)
    >>> plt.plot(x, g.gaussian(x))
    >>> plt.title('Gaussian')
    >>> plt.xlabel('Abcissa')
    >>> plt.ylabel('Ordinate')
    
    >>> # use an array [3] as a single parameter vector
    >>> z = np.array([2., 2, 3])
    >>> plt.plot(x, g.gaussian(x, *z))

    >>> # Test that it integrates to 1.
    >>> a = np.indices([100, 100]) - 50
    >>> print(np.sum(g.gaussian(a, 3, 3)))
    0.999999999999997
    >>> print(np.sum(g.gaussian(a, np.array([1,2]), np.array([2,3]))))
    1.0000000107

    >>> plt.clf()
    >>> plt.imshow(g.gaussian(a, [3,5], [7,3]))
    >>> plt.title('2D Gaussian')
    >>> plt.xlabel('X')
    >>> plt.ylabel('Y')

    >>> A gaussian + a linear background level:
    >>> g2 = g.gaussian(x, width=(1.2, 1.15), center=(13.2,15.75), height=4.3,
    >>>                 bgpars=[0.05, 0.01, 1.0])
    >>> plt.figure(1)
    >>> plt.clf()
    >>> plt.imshow(g2, origin='lower_left', interpolation='nearest')
    >>> plt.colorbar()
    >>> plt.title('2D Gaussian')
    >>> plt.xlabel('X')
    >>> plt.ylabel('Y')

    >>> plt.figure(2)
    >>> plt.clf()
    >>> plt.plot(g2[13,:])
    >>> plt.title('X slice of 2D Gaussian')
    >>> plt.xlabel('X')
    >>> plt.ylabel('Z')

    >>> plt.figure(3)
    >>> plt.clf()
    >>> plt.plot(g2[:,16])
    >>> plt.title('Y slice of 2D Gaussian')
    >>> plt.xlabel('Y')
    >>> plt.ylabel('Z')

    Revisions
    ---------
    2007-09-17 0.1  jh@physics.ucf.edu	Initial version 0.01
    2007-10-02 0.2  jh@physics.ucf.edu	Started making N-dimensional,
    					put width before center in args.
    2007-11-13 0.3  jh@physics.ucf.edu	Fixed docs, bugs, added param,
    					made N-dimensional
    2009-10-01 0.4  jh@physics.ucf.edu	Fixed docs.
    2009-10-25 0.5  jh@physics.ucf.edu	Added examples and plot labels.
    2011-05-03      patricio            Params option no longer sopported,
                                        Added bgpars to add a background.
                                        pcubillos@fulbrightmail.org
  """

  ndim = np.ndim(x) - 1
  if ndim == 0:    # We use an indexing trick below that fails for 1D case.
    ndim = 1
    oldshape = np.shape(x)
    x.shape = (1, x.shape[0])

  # Make center a ndarray:
  if type(center) != np.ndarray:
    center += np.zeros(ndim)

  # Make width a ndarray:
  if type(width) != np.ndarray:
    width  += np.zeros(ndim)
  r2pi = np.sqrt(2. * np.pi)

  # Define height if needed:
  if height == None:
    height = np.product(1. / (width * r2pi))
  ponent = 0.0

  for i in np.arange(ndim):
    ponent += ( (x[i] - center[i]) / width[i] )**2
  if 'oldshape' in locals():
    x.shape = oldshape

  # Set up the background:
  if ndim == 2:
    background = x[0]*bgpars[0] + x[1]*bgpars[1] + bgpars[2]
  else: # it must be 1D:
    background = x*bgpars[0]                     + bgpars[2]

  return height * np.exp(-0.5 * ponent) + background


def old_gaussianguess(y, x=None, mask=None):
  """
    Crudely estimates the parameters of a Gaussian that fits the (y, x) data.

    Parameters
    ----------
    y : ndarray
        Array giving the function values.
    x : ndarray, same shape as y
        (optional) An array of the same shape as y giving the
        abcissas of y (if missing, uses array indices).  Must be
        sorted ascending (which is not checked).

    Returns
    -------
    param : tuple, 3 elements
    This function returns a tuple giving extimates of the (width,
    center, height) of a Gaussian that might fit the input data.
    See 'param' input parameter of gaussian() for format of this
    tuple.

    Notes
    -----
    Currently only works for 1D data.
    
    If the data do not look Gaussian, and certainly if they contain
    spikes higher/lower than the peak of the real Gaussian, the
    parameter estimates will be poor.  x must be sorted ascending
    (which is not checked).

    Method: The most extreme element of y (or its neighbor if a border
    element) is the location and height of the peak.  The routine
    looks for the width at 0.6 of this maximum.

    Todo:
    When expanding to 2D, take arrays of X and Y coords rather than a
    (redundant) 2D array.  This will prevent irregular grids.


    2011-05-05 patricio: This function doesnt work for 2D, I don't
    even know if it works for 1D. If I ever have time I'll see what
    can we do. The function below seems to work fine for our 2D data.

    Examples
    --------
    >>> import gaussian as g

    >>> x = np.arange(-10., 10.05, 0.1)
    >>> y = g.gaussian(x)
    >>> print(g.gaussianguess(y, x))
    (0.99999999999999645, -3.5527136788005009e-14, 0.3989422804014327)


    Revisions
    ---------
    2007-09-17 0.1 jh@physics.ucf.edu    Initial version 0.01
    2007-11-13 0.2 jh@physics.ucf.edu    Fixed docs, return order.
    2008-12-02 0.3 nlust@physics.ucf.edu Fixed a bug where if an
                   initial guess was not provided, it would error out
    2009-10-25 0.4 jh@physics.ucf.edu    Converted to standard doc header.
  """

  if y.ndim != 1 :
    raise ArrayShapeError, "y must be 1D, for now."

  if x == None :
    x = np.indices(y.shape)[0]
  else:
    if x.shape == (1, y.shape):
      oldshape = x.shape
      x.shape  = y.shape
    elif x.shape != y.shape :
      raise ArrayShapeError, "x must have same shape as y (and be sorted)."
      
  # Default mask:
  if mask == None:
    mask = np.ones(np.shape(y))

  ymax  = np.amax(y*mask)
  #iymax = np.where(y == ymax)[0][0]
  iymax = np.argmax(y*mask)

  ymin  = np.amin(y*mask)
  #iymin = np.where(y == ymin)[0][0]
  iymin = np.argmin(y*mask)

  if np.abs(ymin) >= np.abs(ymax):
    icenter = iymin
  else:
    icenter = iymax

  icenter = np.clip(icenter, 1, x.size-2)
  
  center = x[icenter]
  height = y[icenter]
  gtsigma = np.where(y > (0.6 * height))
  width = (x[gtsigma[0].max()] - x[gtsigma[0].min()] ) / 2.

  if 'oldshape' in locals():
    x.shape = oldshape

  return (width, center, height)


def gaussianguess(data, mask=None, yxguess=None):

  # Default mask:
  if mask == None:
    mask = np.ones(np.shape(data))

  # Center position guess, looking the max value:
  if yxguess == None:
    gcenter = np.unravel_index(np.argmax(data*mask), np.shape(data))
  else:
    gcenter = np.around(yxguess[0]), np.around(yxguess[1])

  # Height guess is value at gcenter position:
  gheight = data[gcenter]  

  # The width guess is the sum of the number of pixels that are
  # greater than two sigma of the values in the x and y direction.
  # This gives a (very) rough guess, in pixels, how wide the PSF is.
  sigma = np.array( [np.std(data[:, gcenter[1]]),    # y std (of central column)
                     np.std(data[gcenter[0], :])] )  # x std (of central row)

  gwidth = ( np.sum((data*mask)[:, gcenter[1]] > 2*sigma[0])/2.0,
             np.sum((data*mask)[gcenter[0], :] > 2*sigma[1])/2.0  )

  return (gwidth, gcenter, gheight)


def fitgaussian(y, x=None, bgpars=None, fitbg=0, guess=None,
                mask=None, weights=None, maskg=False, yxguess=None):
  """
    Fits an N-dimensional Gaussian to (value, coordinate) data.

    Parameters
    ----------
    y : ndarray
        Array giving the values of the function.
    x : ndarray
        (optional) Array (any shape) giving the abcissas of y (if
        missing, uses np.indices(y).  The highest dimension must be
        equal to the number of other dimensions (i.e., if x has 6
        dimensions, the highest dimension must have length 5).  The
        rest of the dimensions must have the same shape as y.  Must be
        sorted ascending (which is not checked), if guess is not
        given.
    bgpars : ndarray or tuple, 3-elements
        Background parameters, the elements determine a X- and Y-linearly
        dependant level, of the form:
        f = Y*bgparam[0] + X*bgparam[1] + bgparam[2]
        (Not tested for 1D yet).
    fitbg : Integer
        This flag indicates the level of background fitting:
        fitbg=0: No fitting, estimate the bg as median(data).
        fitbg=1: Fit a constant to the bg (bg = c).
        fitbg=2: Fit a plane as bg (bg = a*x + b*y + c).
    guess : tuple, (width, center, height)
        Tuple giving an initial guess of the Gaussian parameters for
        the optimizer.  If supplied, x and y can be any shape and need
        not be sorted.  See gaussian() for meaning and format of this
        tuple.
    mask : ndarray
        Same shape as y. Values where its corresponding mask value is
        0 are disregarded for the minimization. Only values where the
        mask value is 1 are considered.
    weights : ndarray
        Same shape as y. This array defines weights for the
        minimization, for scientific data the weights should be
        1/sqrt(variance).

    Returns
    -------
    params : ndarray
    
        This array contains the best fitting values parameters: width,
        center, height, and if requested, bgpars. with:
          width :  The fitted Gaussian widths in each dimension.
          center : The fitted Gaussian center coordinate in each dimension.
          height : The fitted height.

    err : ndarray
        An array containing the concatenated uncertainties
        corresponding to the values of params.  For example, 2D input
        gives np.array([widthyerr, widthxerr, centeryerr, centerxerr,
        heighterr]).

    Notes
    -----
    If the input does not look anything like a Gaussian, the result
    might not even be the best fit to that.

    Method: First guess the parameters (if no guess is provided), then
    call a Levenberg-Marquardt optimizer to finish the job.

    Examples
    --------

    >>> import matplotlib.pyplot as plt
    >>> import gaussian as g

    >>> # parameters for X
    >>> lx = -3.  # low end of range
    >>> hx = 5.   # high end of range
    >>> dx = 0.05 # step

    >>> # parameters of the noise
    >>> nc = 0.0  # noice center
    >>> ns = 1.0  # noise width
    >>> na = 0.2  # noise amplitude

    >>> # 1D Example

    >>> # parameters of the underlying Gaussian
    >>> wd = 1.1  # width
    >>> ct = 1.2  # center
    >>> ht = 2.2  # height

    >>> # x and y data to fit
    >>> x  = np.arange(lx, hx + dx / 2., dx)
    >>> x +=                             na * np.random.normal(nc, ns, x.size)
    >>> y  = g.gaussian(x, wd, ct, ht) + na * np.random.normal(nc, ns, x.size)
    >>> s  = x.argsort()    # sort, in case noise violated order
    >>> xs = x[s]
    >>> ys = y[s]

    >>> # calculate guess and fit
    >>> (width, center, height)     = g.gaussianguess(ys, xs)
    >>> (fw,    fc,     fh,    err) = g.fitgaussian(ys, xs)

    >>> # plot results
    >>> plt.clf()
    >>> plt.plot(xs, ys)
    >>> plt.plot(xs,      g.gaussian(xs, wd,    ct,     ht))
    >>> plt.plot(xs,      g.gaussian(xs, width, center, height))
    >>> plt.plot(xs,      g.gaussian(xs, fw,    fc,     fh))
    >>> plt.title('Gaussian Data, Guess, and Fit')
    >>> plt.xlabel('Abcissa')
    >>> plt.ylabel('Ordinate')
    >>> # plot residuals
    >>> plt.clf()
    >>> plt.plot(xs, ys - g.gaussian(xs, fw,    fc,     fh))
    >>> plt.title('Gaussian Fit Residuals')
    >>> plt.xlabel('Abcissa')
    >>> plt.ylabel('Ordinate')

    >>> # 2D Example

    >>> # parameters of the underlying Gaussian
    >>> wd = (1.1, 3.2)  # width
    >>> ct = (1.2, 3.1)  # center
    >>> ht = 2.2         # height

    >>> # x and y data to fit
    >>> nx = (hx - lx) / dx + 1
    >>> x  = np.indices((nx, nx)) * dx + lx
    >>> y  = g.gaussian(x, wd, ct, ht) + na * np.random.normal(nc, ns, x.shape[1:])

    >>> # calculate guess and fit
    >>> #(width, center, height)     = g.gaussianguess(y, x) # not in 2D yet...
    >>> (fw,    fc,     fh,    err) = g.fitgaussian(y, x, (wd, ct, ht))

    >>> # plot results
    >>> plt.clf()
    >>> plt.title('2D Gaussian Given')
    >>> plt.xlabel('X')
    >>> plt.ylabel('Y')
    >>> plt.imshow(    g.gaussian(x, wd,    ct,     ht))
    >>> plt.clf()
    >>> plt.title('2D Gaussian With Noise')
    >>> plt.xlabel('X')
    >>> plt.ylabel('Y')
    >>> plt.imshow(y)
    >>> #plt.imshow(    g.gaussian(x, width, center, height)) # not in 2D yet...
    >>> plt.clf()
    >>> plt.title('2D Gaussian Fit')
    >>> plt.xlabel('X')
    >>> plt.ylabel('Y')
    >>> plt.imshow(    g.gaussian(x, fw,    fc,     fh))
    >>> plt.clf()
    >>> plt.title('2D Gaussian Fit Residuals')
    >>> plt.xlabel('X')
    >>> plt.ylabel('Y')
    >>> plt.imshow(y - g.gaussian(x, fw,    fc,     fh))

    >>> # All cases benefit from...

    >>> # show difference between fit and underlying Gaussian
    >>> # Random data, your answers WILL VARY.
    >>> np.array(fw) - np.array(wd)
    array([ 0.00210398, -0.00937687])
    >>> np.array(fc) - np.array(ct)
    array([-0.00260803,  0.00555011])
    >>> np.array(fh) - np.array(ht)
    0.0030143371034774269

    >>> Last Example:
    >>> x = np.indices((30,30))
    >>> g1 = g.gaussian(x, width=(1.2, 1.15), center=(13.2,15.75), height=1e4,
    >>>                 bgpars=[0.0, 0.0, 100.0])
    >>> error = np.sqrt(g1) * np.random.randn(30,30)
    >>> y = g1 + error
    >>> var = g1
    >>> 
    >>> plt.figure(1)
    >>> plt.clf()
    >>> plt.imshow(y, origin='lower_left', interpolation='nearest')
    >>> plt.colorbar()
    >>> plt.title('2D Gaussian')
    >>> plt.xlabel('X')
    >>> plt.ylabel('Y')
    >>> 
    >>> guess = ((1.2,1.2),(13,16.),1e4)
    >>> reload(g)
    >>> fit = g.fitgaussian(y, x, bgpars=[0.0, 0.0, 110.], fitbg=1, guess=guess,
    >>>                     mask=None, weights=1/np.sqrt(var))
    >>> print(fit[0])


    Revisions
    ---------
    2007-09-17  Joe      Initial version, portions adapted from
                         http://www.scipy.org/Cookbook/FittingData.
                         jh@physics.ucf.edu
    2007-11-13  Joe      Made N-dimensional.
    2008-12-02  Nate     Included error calculation, and return Fixed a bug
                         in which if the initial guess was None, and incorrect
                         shape array was generated. This caused gaussian guess
                         to fail.
                         nlust@physics.ucf.edu
    2009-10-25           Converted to standard doc header, fixed examples to
                         return 4 parameters.
    2011-05-03  patricio Added mask, weights, and background-fitting options.
                         pcubillos@fulbrightmail.org
  """

  if x == None:
    x = np.indices(np.shape(y))
  else:
    if (   ((x.ndim == 1) and (x.shape     != y.shape))
        or ((x.ndim >  1) and (x.shape[1:] != y.shape))):
      raise ValueError, "x must give coordinates of points in y."

  # Default mask: all good
  if mask == None:
    mask    = np.ones(np.shape(y))

  # Default weights: no weighting
  if weights == None:
    weights = np.ones(np.shape(y))

  # Mask the gaussian if requested:
  medmask = np.copy(mask)
  if maskg and (yxguess != None or guess != None):
    if yxguess != None:
      center = yxguess
    elif guess != None:
      center = guess[1]
    medmask *= (1 - d.disk(3, center, np.shape(y)))

  # Estimate the median of the image:
  medbg = np.median(y[np.where(medmask)])

  if bgpars == None:
    bgpars = [0.0, 0.0, medbg]

  # get a guess if not provided
  if guess == None:
    if yxguess == None:
      guess = gaussianguess(y-medbg, mask=mask)
    else:
      guess = gaussianguess(y-medbg, mask=mask, yxguess=yxguess)

  # "ravel" the guess
  gparams = np.append(guess[0], guess[1])
  gparams = np.append(gparams,  guess[2])

  # Background params to fit:
  if fitbg == 0:
    bgparams = []
  elif fitbg == 1:
    bgparams = bgpars[2]
  elif fitbg == 2:
    bgparams = bgpars

  # Concatenate sets of parameters we want to fit:
  params = np.append(gparams, bgparams)
  # Rest of parameters needed by residuals:
  args = (x, y, mask, weights, bgpars, fitbg)
  
  # The fit:
  p, cov, info, mesg, success = so.leastsq(residuals, params, args,
                                           full_output=True)
  try:
    err = np.sqrt(np.diagonal(cov))
  except:
    err = None
  
  return p, err


def residuals(params, x, data, mask, weights, bgpars, fitbg):
  """
    Calculates the residuals between data and a gaussian model
    determined by the rest of the parameters. Used in fitgaussian.
    
    Parameters
    ----------
    params : 1D ndarray
        This array contains the parameters desired to fit with
        fitgaussian, depending on fitbg, the number of elements
        varies.
    x : ndarray
        Array (any shape) giving the abcissas of data.
    data : ndarray
        Array giving the values of the function.
    mask : ndarray
        Same shape as data. Values where its corresponding mask value is
        0 are disregarded for the minimization. Only values where the
        mask value is 1 are considered.
    weights : ndarray
        Same shape as data. This array defines weights for the
        minimization, for scientific data the weights should be
        1/sqrt(variance).
    bgpars : ndarray or tuple, 3-elements
        Background parameters, the elements determine a X- and Y-linearly
        dependant level, of the form:
        background = Y*bgparam[0] + X*bgparam[1] + bgparam[2]
    fitbg : Integer
        This flag indicates the level of background fitting:
        fitbg=0: No fitting, estimate the bg as median(data).
        fitbg=1: Fit a constant to the bg (bg = c).
        fitbg=2: Fit a plane as bg (bg = a*x + b*y + c).

    Returns
    -------
    residuals : 1D ndarray
        An array of the (unmasked) weighted residuals between data and
        a gaussian model determined by params (and bgpars when
        necessary).

    Examples
    --------
    
    Revisions
    ---------
    2011-05-03  patricio Initial version.
                         pcubillos@fulbrightmail.org
  """
  
  # Use bgpars as default for background parameters, if those values
  # are being fitted update them:
  bgparams = bgpars
  if fitbg == 1:
    bgparams[2] = params[-1]   # update
    params      = params[0:-1] # remove last parameters from params
  elif fitbg == 2:
    bgparams = params[-3:]   # update
    params   = params[0:-3]  # remove last parameters

  # Extract width, center, and height from params:
  data_dims = np.ndim(data)
  params_len = len(params)

  width  = params[0        :  data_dims]
  center = params[data_dims:2*data_dims]
  if params_len - 2*data_dims == 1:
    height = params[2*data_dims]
  else:
    # when height is None, queda la cagada, avoid this case.
    height = None

  # Produce the model:
  model = gaussian(x, width, center, height, bgparams).squeeze()
  # Calculate residuals:
  res = (model - data) * weights
  # Return only unmasked values:
  return res[np.where(mask)]
  
  
def gaussians(x, param):
  """
    Evaluate more than 1 gaussian.
  """

  ndim = x.ndim - 1
  if ndim == 0:    # We use an indexing trick below that fails for 1D case.
    ndim = 1
    oldshape = x.shape
    x.shape = (1, x.shape[0])

  # The number of gaussians:
  ngauss = np.shape(param)[0]
  if ngauss == 1:
    param = [param]

  result = np.zeros(x[0].shape)
  for k in np.arange(ngauss):  # Unpack parameters
    pdim = len(param[k])
    if pdim % 2:  # pdim is odd (when height is specified)
      pdim = (pdim - 1) / 2
      height = param[k][-1]
    else:         # pdim is even
      pdim = pdim / 2
      height = None
    width  = param[k][     :     pdim]
    center = param[k][pdim : 2 * pdim]

    if type(center) != np.ndarray:
      center += np.zeros(ndim)
    if type(width) != np.ndarray:
      width  += np.zeros(ndim)
    if height == None:
      height = np.product(1.0 / (width * np.sqrt(2.0 * np.pi)))
    ponent = 0.0
    for i in np.arange(ndim):
      ponent += ( (x[i] - center[i]) / width[i] )**2.0
    result += height * np.exp(-0.5 * ponent)

  if 'oldshape' in locals():  # reshape it back if necessary 
    x.shape = oldshape
  return result


def fitgaussians(y, x=None, guess=None, sigma=1.0):
  """
    Fit position and flux of a data image with gaussians, same sigma
    is applied to all dispersions.


    Parameters:
    -----------
    y : array_like
        Array giving the values of the function.
    x : array_like
        (optional) Array (any shape) giving the abcissas of y (if
        missing, uses np.indices(y).
        
    guess : 2D-tuple, [[width1, center1, height1],
                       [width2, center2, height2],
                       ...                       ]
        Tuple giving an initial guess of the Gaussian parameters for
        the optimizer.  If supplied, x and y can be any shape and need
        not be sorted.  See gaussian() for meaning and format of this
        tuple.

  """
  if x == None:
    x = np.indices(y.shape)[0]
  else:
    if (   ((x.ndim == 1) and (x.shape     != y.shape))
        or ((x.ndim >  1) and (x.shape[1:] != y.shape))):
      raise ValueError, "x must give coordinates of points in y."

  # "ravel" the guess
  ngauss = np.shape(guess)[0]
  params = np.ravel(guess)
  params = np.append(guess, sigma)

  # Minimize residuals of the fit:
  p, cov, info, mesg, success = so.leastsq(resids, params, args=(x,ngauss,y),
                                           full_output=True)

  sigma = p[-1]
  p   = np.reshape(p  [0:-1], (ngauss, len(p  [0:-1])/ngauss))

  iscov = 0 if cov==None else 1
  extra = (p, sigma, iscov, cov, info, mesg)

  return np.array(p[0,0:2]), extra


def resids(param, x, ngauss, y):
  sigma = param[-1]
  param = np.reshape(param[0:-1], (ngauss, len(param[0:-1])/ngauss))

  gss = []
  for k in np.arange(ngauss):
    gauss = np.append(sigma, np.append(sigma, param[k]))
    gss = np.append(gss,gauss)
  p = np.reshape(gss, (ngauss,len(gss)/ngauss))
  return np.ravel(gaussians(x,param=p)-y)

