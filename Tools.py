import astropy.io.fits

def readFitsData(filename, verbose=False):
  '''Read in data from a FITS image.'''
  hdu = astropy.io.fits.open(filename)
  if verbose:
    print "      read ", filename
  image = hdu[0].data
  hdu.close()
  return image

def writeFitsData(data, filename, verbose=False):
  '''Write data to a FITS image.'''
  hdu = astropy.io.fits.PrimaryHDU(data)
  hdu.writeto(filename, clobber=True)
  if verbose:
    print "      wrote image to ", filename
