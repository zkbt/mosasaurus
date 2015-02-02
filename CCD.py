from imports import *

class CCD(Talker):
  '''CCD object handles every related to an individual CCD exposure.'''

  def __init__(self, obs, n=0, type=None, calib=None, verbose=True, **kwargs):
    '''Initialized the CCD object.'''


    # decide whether or not this Reducer is chatty
    Talker.__init__(self, **kwargs)

    # initialize the basic stuff we need
    self.obs = obs
    self.calib = calib
    self.display = zachopy.display.ds9('CCD')
    self.space = '  '



    # set up options
    self.flags = {'subtractbias':True, 'subtractdark':True, 'multiplygain':False}
    self.visualize = True
    self.verbose = verbose

    # point this CCD to right file, if provided
    self.set(n, type)

    # make sure the stitched directory is defined
    zachopy.utils.mkdir(self.obs.workingDirectory+'/stitched/'	)


  def set(self, n=None, type=None):
    '''Point this CCD object at a specific file.'''

    # keep track if this is some special kind of image
    self.type = type

    # define the file prefix
    if n is not None:
      self.fileprefix = self.obs.fileprefix(n)

      # define a nickname for referring to this image
      if self.type is None:
        label = 'ccd'
      else:
        label = self.type
      self.name = '{1}{0:04d}'.format(n,label)

      # define a stitched filename
      self.stitched_filename = self.obs.workingDirectory + 'stitched/{0}.fits'.format(self.name)

    # empty out the header and data variables
    self.header = None
    self.data = None

    # print status
    if self.verbose:
      self.speak(self.space + "set image to {0}".format(self.name))

  def readHeader(self, n=None, type=None):
    '''Read in the header for this image.'''

    # set the CCD, if necessary
    if n is not None:
      self.set(n, type)

    # read one of the two images to get header
    filename = self.fileprefix + 'c1.fits'
    hdu = astropy.io.fits.open(filename)
    header = hdu[0].header
    hdu.close()

    # store the header in this object
    self.header = header

    if self.verbose:
      self.speak(self.space + "read header from {0}".format(self.name))

    # return the header
    return self.header

  def readData(self, n=None, type=None):
    '''Read in the image data for this exposure (create stitched exposure if necessary.)'''

    # set the CCD, if necessary
    if n is not None:
      self.set(n, type)

    # try loading a stitched image, otherwise create one from scratch
    try:
      self.data = readFitsData(self.stitched_filename)
    except:
      self.createStitched()

    # print status
    if self.verbose:
      self.speak(self.space + "read image from {0}".format(self.name))

    # return the image
    return self.data

  def loadOverscanTrimHalfCCD(self, filename):
    '''Open one half of an LDSS3 CCD, subtract the overscan, and trim.'''

    # open the FITS file, split into header and data
    hdu = astropy.io.fits.open(filename)
    header = hdu[0].header
    data = readFitsData(filename)

    # take the parts of CCD exposed to light
    goodData = data[self.obs.databottom:self.obs.datatop,self.obs.dataleft:self.obs.dataright]
    goodBias = data[self.obs.databottom:self.obs.datatop,self.obs.dataright:]

    # estimate the 1D bias (and drawdown, etc...) from the overscan
    biasEstimate = np.median(goodBias, axis=1)
    biasImage = np.ones(goodData.shape)*biasEstimate[:,np.newaxis]

    return (goodData - biasImage), header

  def createStitched(self):
    '''Create and load a stitched CCD image, given a file prefix.'''

    # print status
    if self.verbose:
      self.speak(self.space + "creating a stitched image for {0}".format( self.stitched_filename))

    # provide different options for different kinds of images
    if self.type == 'Bias':
      self.flags['subtractbias'] = False
      self.flags['subtractdark'] = False
      self.flags['multiplygain'] = False
    elif self.type == 'Dark':
      self.flags['subtractbias'] = True
      self.flags['subtractdark'] = False
      self.flags['multiplygain'] = False
    elif self.type == 'FlatInADU':
      self.flags['subtractbias'] = True
      self.flags['subtractdark'] = True
      self.flags['multiplygain'] = False
    else:
      self.flags['subtractbias'] = True
      self.flags['subtractdark'] = True
      self.flags['multiplygain'] = True

    # don't restitch if unnecessary
    if os.path.exists(self.stitched_filename) and self.verbose:
      self.speak(self.space + "{0} has already been stitched".format(self.name))
    else:
      # process the two halves separately, and then smush them together
      filenames = [self.fileprefix + 'c1.fits', self.fileprefix + 'c2.fits']

      # load the two halves
      c1data, c1header = self.loadOverscanTrimHalfCCD(filenames[0])
      c2data, c2header = self.loadOverscanTrimHalfCCD(filenames[1])

      if self.visualize:
        tempstitched = np.hstack([c1data, np.fliplr(c2data)])

      # divide by the gain (KLUDGE! make sure these are the best estimates!)
      if self.flags['multiplygain']:

        try:
          self.calib.gains
        except:
          self.calib.estimateGain()
        self.speak("       multiplying by gains of {0} e-/ADU".format(self.calib.gains))
        c1data *= self.calib.gains[0]
        c2data *= self.calib.gains[1]

      # stitch the CCD's together
      stitched = np.hstack([c1data, np.fliplr(c2data)])
      self.speak("       stitching images of size {0} and {1} into one {2} image".format(c1data.shape, c2data.shape, stitched.shape))

      # subtract bias
      if self.flags['subtractbias']:
        self.speak("       subtracting bias image")
        stitched -= self.calib.bias()

      # normalize darks by exposure time
      if self.type == 'Dark':
        stitched /= c1header['EXPTIME']

      # subtract dark
      if self.flags['subtractdark']:
        self.speak("       subtracting dark image")
        stitched -= self.calib.dark()*c1header['EXPTIME']



      # save the stitched image into memory
      self.data = stitched

      # write out the image to a stitched image
      writeFitsData(stitched, self.stitched_filename)
      if self.verbose:
        self.speak(self.space + "stitched and saved {0}".format(self.name))

      # if need be, visualize the images
      if False:#self.visualize:
        self.display.many([stitched, tempstitched])

  def amplifiers(self):
    return self.data[:,0:self.obs.dataright - self.obs.dataleft], self.data[:,self.obs.dataright - self.obs.dataleft:]


  def loadImages(self, n, type=None):
    '''Load a series of CCD images, returning them as a cube.'''

    # if n is a single element array, just return one image
    try:
      n[1]
    except:
      self.set(n[0])
      return self.readData(type=type)

    # loop over the image numbers, and read them
    images = []
    for i in range(len(n)):
      images.append(self.readData(n[i], type))

    # convert a list of images into an array
    cube = np.array(images)

    # return that array
    return cube
