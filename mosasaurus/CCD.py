from imports import *

class CCD(Talker):
  '''CCD object handles every related to an individual CCD exposure.'''

  def __init__(self, obs, n=0, imageType=None, calib=None, verbose=True, **kwargs):
    '''Initialized the CCD object.'''


    # decide whether or not this Reducer is chatty
    Talker.__init__(self, **kwargs)

    # initialize the basic stuff we need
    self.obs = obs
    self.calib = calib
    self.display = self.calib.display

    self.space = '  '



    # set up options
    self.flags = {'subtractbias':True, 'subtractdark':True, 'multiplygain':False}
    self.visualize = False
    self.verbose = verbose

    # point this CCD to right file, if provided
    self.set(n, imageType)

    # make sure the stitched directory is defined
    zachopy.utils.mkdir(self.obs.workingDirectory+'/stitched/'	)


  def set(self, n=None, imageType=None):
    '''Point this CCD object at a specific file.'''

    # keep track if this is some special kind of image
    self.imageType = imageType
    self.n = n

    # define the file prefix
    if n is not None:
      self.fileprefix = self.obs.fileprefix(n)

      # define a nickname for referring to this image
      if self.imageType is None:
        label = 'ccd'
      else:
        label = self.imageType
      self.name = '{1}{0:04d}'.format(n,label)

      # define a stitched filename
      self.stitched_filename = self.obs.workingDirectory + 'stitched/{0}.fits'.format(self.name)

    # empty out the header and data variables
    self.header = None
    self.data = None

    # print status
    #if self.verbose:
    #  self.speak(self.space + "set image to {0}".format(self.name))

  def readHeader(self, n=None, imageType=None):
    '''Read in the header for this image.'''

    # set the CCD, if necessary
    if n is not None:
      self.set(n, imageType)

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

  def writeData(self):
      self.speak('saving image to {0}'.format(self.stitched_filename))

      # debugging!
      #if self.imageType == 'Science':
      #      self.display.one(self.data)
      #      self.input('cosmics or no?')
      writeFitsData(self.data, self.stitched_filename)

  def readData(self, n=None, imageType=None):
    '''Read in the image data for this exposure (create stitched exposure if necessary.)'''

    # set the CCD, if necessary
    if n is not None:
      self.set(n, imageType)

    # try loading a stitched image, otherwise create one from scratch
    try:
      self.data = readFitsData(self.stitched_filename)
    except IOError:
      self.createStitched()

    if imageType == 'Science':
        self.cosmicdiagnostic = np.load(self.obs.workingDirectory+'cosmics/rejectedpercolumn{0:04.0f}.npy'.format(n))

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
    if self.imageType == 'Bias':
      self.flags['subtractbias'] = False
      self.flags['subtractdark'] = False
      self.flags['multiplygain'] = False
      self.flags['subtractcrosstalk'] = False
    elif self.imageType == 'Dark':
      self.flags['subtractbias'] = True
      self.flags['subtractdark'] = False
      self.flags['multiplygain'] = False
      self.flags['subtractcrosstalk'] = False
    elif self.imageType == 'FlatInADU':
      self.flags['subtractbias'] = True
      self.flags['subtractdark'] = True
      self.flags['multiplygain'] = False
      self.flags['subtractcrosstalk'] = False
    else:
      self.flags['subtractbias'] = True
      self.flags['subtractdark'] = True
      self.flags['multiplygain'] = True
      self.flags['subtractcrosstalk'] = True

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



      if self.flags['subtractcrosstalk']:
          # is this possible?
          pass

      # stitch the CCD's together
      stitched = np.hstack([c1data, np.fliplr(c2data)])
      #self.speak("stitching images of size {0} and {1} into one {2} image".format(c1data.shape, c2data.shape, stitched.shape))

      if self.visualize:
          self.display.one(stitched, clobber=True)
          self.input('This is the raw stitched image; press enter to continue.')

      # subtract bias
      if self.flags['subtractbias']:
        self.speak("subtracting bias image")
        stitched -= self.calib.bias()

      if self.visualize:
          self.display.one(stitched, clobber=True)
          self.input('after subtracting bias')

      # normalize darks by exposure time
      if self.imageType == 'Dark':
        stitched /= c1header['EXPTIME']

      # subtract dark
      if self.flags['subtractdark']:
        self.speak("subtracting dark image")
        stitched -= self.calib.dark()*c1header['EXPTIME']

      if self.visualize:
          self.display.one(stitched, clobber=True)
          self.visualize = self.input('after subtracting dark; type [s] to stop showing these').lower() != 's'

      # divide by the gain (KLUDGE! make sure these are the best estimates!)
      if self.flags['multiplygain']:

        try:
          self.calib.gains
        except:
          self.calib.estimateGain()
        self.speak("multiplying by gains of {0} e-/ADU".format(self.calib.gains))
        gain1 = np.zeros_like(c1data) + self.calib.gains[0]
        gain2 = np.zeros_like(c2data)+ self.calib.gains[1]
        gainimage = np.hstack([gain1, np.fliplr(gain2)])
        stitched *= gainimage

      if self.visualize:
          self.display.one(stitched, clobber=True)
          self.visualize = self.input('after multiplying by gain; type [s] to stop showing these').lower() != 's'

      # put the stitched image into the CCD's memory
      self.data = stitched

      if self.imageType == 'Science':
          self.rejectCosmicRays()

      # write out the image to a stitched image
      writeFitsData(self.data, self.stitched_filename)
      if self.verbose:
        self.speak(self.space + "stitched and saved {0}".format(self.name))

  def rejectCosmicRays(self, remake=False, threshold=7.5, visualize=False, nBeforeAfter=5):
     '''Stitch all science images, establish a comparison noise level for each pixel.'''
     # make sure a cosmics directory exists
     cosmics_directory = self.obs.workingDirectory + 'cosmics/'
     zachopy.utils.mkdir(cosmics_directory)

     # figure out how many images to consider
     nImages = 2*nBeforeAfter + 1
     imageType = 'ScienceUnmitigated'
     n = self.n
     try:
         #print "testing"
         ok = remake == False
         #print 'remake'
         ok = ok & os.path.exists(self.obs.workingDirectory + 'stitched/Science{0:04.0f}.fits'.format(n))
         #print 'fits'
         ok = ok & os.path.exists(self.obs.workingDirectory + 'cosmics/rejectedpercolumn{0:04.0f}.npy'.format(n))
         #print 'rejected'
         #print ok
         assert(ok)
         self.speak('a cosmic-rejected stitched/Science{0:04.0f}.fits already exists!'.format(n))
     except (IOError,AssertionError):
         nComparison = np.arange(-nBeforeAfter + n, nBeforeAfter+n+1, 1)
         nComparison = nComparison[(nComparison >= np.min(self.obs.nScience)) & (nComparison <= np.max(self.obs.nScience))]
         comparison = self.loadImages(n=nComparison, imageType=imageType)
         shape = np.array(comparison.shape)
         axis =0
         shape[axis] = 1

         self.speak('comparing image {0} to images {1} to remove cosmic rays'.format(n, nComparison))
         image = self.data
         #image = comparison[nComparison == n,:,:].squeeze()


         # calculate median and noise of comparisons
         med = np.median(comparison, axis=axis)
         noise = np.maximum(1.48*np.median(np.abs(comparison - med.reshape(shape)), axis=axis), 1.0)

         bad = (image - med)/noise > threshold
         corrected = image + 0.0
         corrected[bad] = med[bad]
         if visualize:
             self.display.replace(image,0)
             self.display.replace(image - corrected,1)
             self.display.replace(corrected,2)
             self.display.scale(mode='zscale')
             self.display.match()
             self.display.tile('column')


         images = [corrected]
         labels = ['Science']
         for i in np.arange(len(labels)):
             self.set(n, labels[i])
             self.data = images[i]
             # self.writeData()

         self.speak('total corrected flux is {0}'.format(np.sum(image - corrected)))

         lostflux = np.sum(image - corrected, axis=0)
         try:
             self.cosmicplot.set_ydata(lostflux)
         except AttributeError:
             plt.figure('cosmic ray rejection', figsize=(5, 3), dpi=100)
             self.axcr = plt.subplot()
             self.axcr.cla()
             self.cosmicplot = self.axcr.plot(lostflux, color='Sienna')[0]
             self.axcr.set_ylim(1.0, 1e8)
             self.axcr.set_xlim(-1, len(lostflux)+1)
             self.axcr.set_yscale('log')
             self.axcr.set_ylabel('Total Flux Rejected (e-/column)')
             self.axcr.set_xlabel('Column (pixels)')
             plt.tight_layout()

         self.axcr.set_title('Science{0:04.0f}'.format(n))
         plt.draw()
         np.save(cosmics_directory + 'rejectedpercolumn{0:04.0f}.npy'.format(n), lostflux)
         plt.savefig(cosmics_directory + 'rejectedpercolumn{0:04.0f}.png'.format(n))
         self.speak('saved cosmic ray rejection checks to {0}'.format(cosmics_directory))
         #self.input('thoughts on CR?')



  def amplifiers(self):
    return self.data[:,0:self.obs.dataright - self.obs.dataleft], self.data[:,self.obs.dataright - self.obs.dataleft:]


  def loadImages(self, n, imageType=None):
    '''Load a series of CCD images, returning them as a cube.'''

    # if n is a single element array, just return one image
    try:
      n[1]
    except:
      self.set(n[0])
      return self.readData(imageType=imageType)

    # loop over the image numbers, and read them
    images = []
    for i in range(len(n)):
      images.append(self.readData(n[i], imageType))

    # convert a list of images into an array
    cube = np.array(images)

    # return that array
    return cube
