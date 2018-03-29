# these are some general tools, to be used in multiple parts of mosasaurus
import astropy.io.fits, os, numpy as np

def readFitsData(filename, verbose=False):
    '''Read in data from a FITS image (ignoring the header).'''
    hdu = astropy.io.fits.open(filename)
    if verbose:
        print("      read ", filename)
    image = hdu[0].data
    hdu.close()
    return image

def writeFitsData(data, filename, verbose=False):
    '''Write data to a FITS image (ignoring the header).'''
    hdu = astropy.io.fits.PrimaryHDU(data.astype(np.float32))
    hdu.writeto(filename, clobber=True)
    if verbose:
        print("      wrote image to ", filename)

def truncate(str, n=12, etc=' ...'):
	'''If a string is too long, truncate it with an "etc..."'''
	if len(str) > n:
		return str[0:n-len(etc)] + etc
	else:
		return ("{0: ^%d}" % n).format(str)

def mkdir(path):
	'''A mkdir that doesn't complain if it fails.'''
	try:
		os.mkdir(path)
	except:
		pass

def mad(x):
	'''Median absolute deviation from the median.'''
	med = np.median(x)
	return np.median(np.abs(x - med))

# keyword arguments for writing astropy tables
tablekw = dict(format='ascii.fixed_width', delimiter='|', bookend=False)
