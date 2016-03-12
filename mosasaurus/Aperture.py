from imports import *
from Tools import *
from Trace import Trace
from WavelengthCalibrator import WavelengthCalibrator

plt.ion()
class Aperture(Talker):
  '''Aperture objects handle individual apertures (e.g. slits),
        storing data and handling extraction. '''

  def __init__(self, x, y, mask, **kwargs):
    '''Initialize an aperture, taking a position on the detector as input.'''

    # decide whether or not this Reducer is chatty
    Talker.__init__(self, **kwargs)
    self.visualize = True
    self.mask = mask
    self.calib = self.mask.calib
    self.obs = self.calib.obs
    self.display = self.calib.display
    #zachods9('aperture', xsize=800, ysize=200, rotate=90)
    self.setup(x,y)
    self.createCalibStamps()
    # testing!
    #self.trace.fit()
    self.createTrace()
    self.createWavelengthCal()

  def stampFilename(n):
	'''Spit out the right stamp filename for this aperture, cut out of CCDn.'''
	return self.directory + 'stamp{0:04}.fits'.format(n)

  def setup(self,x,y):
    '''Setup the basic geometry of the aperture.'''
    if self.obs.instrument == 'LDSS3C':
      self.x = x
      self.y = y
      self.maskWidth = self.obs.subarray
      #(self.obs.skyWidth + self.obs.skyGap)*2 #+20 +
      self.ystart = np.maximum(y - self.obs.blueward, 0)
      self.yend = np.minimum(y + self.obs.redward, self.obs.ysize)
      self.xstart = np.maximum(x - self.maskWidth, 0)
      self.xend = np.minimum(x + self.maskWidth, self.obs.xsize)
      # remember python indexes arrays by [row,column], which is opposite [x,y]
      x_fullframe, y_fullframe = np.meshgrid(np.arange(self.calib.images['Science'].shape[1]),
                        np.arange(self.calib.images['Science'].shape[0]))
      self.x_sub = x_fullframe[self.ystart:self.yend, self.xstart:self.xend]
      self.y_sub = y_fullframe[self.ystart:self.yend, self.xstart:self.xend]

      self.xbox = (self.xstart + self.xend)/2
      self.ybox = (self.ystart + self.yend)/2
      self.wbox = np.abs(self.xend - self.xstart)
      self.hbox = np.abs(self.yend - self.ystart)

      # first index of np. array is in the wavelength (w) direction
      # second index is in the spatial (s) direction
      # we'll define these now to help keep things straight
      self.w = self.y_sub - self.y#self.ystart
      self.s = self.x_sub - self.x#self.xstart
      self.windex = 0
      self.sindex = 1 - self.windex
      if self.windex == 0:
        self.waxis = self.w[:,0]
        self.saxis = self.s[0,:]


      self.name = 'aperture_{0:.0f}_{1:.0f}'.format(self.x, self.y)
      self.directory = self.obs.extractionDirectory + self.name + '/'
      zachopy.utils.mkdir(self.directory)
      self.speak("created a spectroscopic aperture at ({0:.1f}, {1:.1f})".format(self.x, self.y))
    else:
      self.speak("*!$!%()@! no cameras besides LDSS3C have been defined yet!")


  def stamp(self, image):
    '''Return a postage stamp of an image, appropriate for this aperture.'''
    return image[self.ystart:self.yend, self.xstart:self.xend]

  def createCalibStamps(self, visualize=True, interactive=True):
    '''Populate the necessary postage stamps for the calibrations.'''
    self.speak("populating calibration stamps")
    filename = self.directory + 'calibStamps_{0}.npy'.format(self.name)
    try:
      self.images = np.load(filename)[()] # i don't understand why i need the "empty tuple" but the internet says so
      self.speak("loaded calibration stamps from {0}".format(filename))
    except:
      self.speak("cutting them for the first time out of the stitched master images")

      # define an empty dictionary to store the calibration stamps
      self.images = {}

      # include the coordinate system over the grid of the image
      self.images['s'] = self.s
      self.images['w'] = self.w
      self.speak("populating calibration stamps with the (s)patial and (w)avelength pixel coordinate images")


      # cut out stamps from the big images
      interesting = ['Science' ,
                    'WideFlat',
                    'He', 'Ne', 'Ar',
                    'BadPixels', 'Dark', 'Bias']
      for k in interesting:
        self.images[k] = self.stamp(self.calib.images[k])

      #self.speak('these are the stamps before interpolating over bad pixels')
      #self.displayStamps(self.images, keys = ['Science', 'WideFlat', 'BadPixels'])
      #self.input('', prompt='(press return to continue)')
      for k in interesting:
          if k != 'BadPixels':
              self.images[k] = zachopy.twod.interpolateOverBadPixels(self.images[k], self.images['BadPixels'])

      #self.speak('and these are they after interpolating over bad pixels')
      #self.displayStamps(self.images, keys = ['Science', 'WideFlat', 'BadPixels'])
      #self.input('', prompt='(press return to continue)')

      # subtract dark from everything but the dark
      #for k in self.images.keys():
      #	if k is not 'Dark':
      #		self.images[k] -= self.images['Dark']

      # create a normalized flat field stamp, dividing out the blaze + spectrum of quartz lamp
      raw_flatfield = self.images['WideFlat']
      overbig_flatfield = np.ones_like(raw_flatfield)
      envelope = np.median(raw_flatfield, self.sindex)
      n_envp = 30
      points = np.linspace(np.min(self.waxis), np.max(self.waxis),n_envp+2)
      spline = scipy.interpolate.LSQUnivariateSpline(self.waxis,envelope,points[1:-2],k=2)
      self.images['NormalizedFlat'] = self.images['WideFlat']/spline(self.waxis).reshape((self.waxis.shape[0],1))
      self.images['NormalizedFlat'] /= np.median(self.images['NormalizedFlat'], self.sindex).reshape(self.waxis.shape[0], 1)




      '''
      if self.visualize:
        try:
            self.ax.cla()
        except:
            self.fi, self.ax = plt.subplots(1,1)
        self.ax.cla()
        self.ax.plot(self.waxis, envelope)
        self.ax.plot(self.waxis, spline(self.waxis))

      if self.visualize:
        self.displayStamps(self.images, keys = ['Science', 'Sky', 'Subtracted', 'WideFlat', 'NormalizedFlat'])
        assert("n" not in self.input('Do you like the calibration stamps?').lower())

      '''
      np.save(filename, self.images)
      self.speak("saved calibration stamps to {0}".format( filename))

  def displayStamps(self, images, keys=None):
    '''Display stamps relevant to this aperture in ds9.'''
    if keys is None:
        keys = images.keys()

    self.display.tile('column')
    for i in range(len(keys)):
      self.display.replace(images[keys[i]], i)

  def createTrace(self):
    '''Fit for the position and width of the trace.'''

    self.speak("populating the trace parameters")
    tracefilename = self.directory + 'trace_{0}.npy'.format(self.name)
    skyfilename = self.directory + 'skyMask_{0}.npy'.format(self.name)
    try:
        traceCoeff, width = np.load(tracefilename)
        self.speak("loaded trace from {0}".format(tracefilename))

        self.images['skyMask'] = np.load(skyfilename)
        self.speak("loaded sky apertures from {0}".format(skyfilename))
    except IOError:
        self.trace = Trace(self)
        traceCoeff, width = self.trace.tracefitcoef, self.trace.tracefitwidth
        self.images['skyMask'] = self.trace.skymask

    self.traceCenter = np.poly1d(traceCoeff)
    self.traceWidth = width
    self.images['RoughLSF'] = np.exp(-0.5*((self.s - self.traceCenter(self.w))/self.traceWidth)**2)
    self.displayTrace()

  def displayTrace(self):
      self.display.rgb( self.images['Science'],
                        self.images['RoughLSF'],
                        self.images['skyMask'])


  def createWavelengthCal(self, remake=False):
    '''Populate the wavelength calibration for this aperture.'''
    self.speak("populating wavelength calibration")
    self.wavelengthcalibrator = WavelengthCalibrator(self)

  def ones(self):
    '''Create a blank array of ones to fill this aperture.'''
    return np.ones_like(self.images['Science'])

  @property
  def extractedFilename(self):
      return self.directory + 'extracted{0:04}.npy'.format(self.n)

  def extract(self, n=0, image=None, subtractsky=True, arc=False, cosmics=False, remake=False):
    '''Extract the spectrum from this aperture.'''
    self.n = n
    try:
        assert(remake == False)
        assert(arc == False)
        self.extracted = np.load(self.extractedFilename)
        self.speak('loaded extracted spectrum from {0}'.format(self.extractedFilename))
    except:

        # make sure we don't subtract the sky, if the arcs are set
        subtractsky = subtractsky & (arc == False)

        # create a dictionary to store the self.extracted spectra
        self.extracted = {}
        self.extracted['w'] = self.waxis

        if image is None:
            # make sure the right data have been loaded
            assert(self.mask.ccd.n == n)

            # extract the raw science stamp from this mask image
            raw = self.stamp(self.mask.ccd.data)
        else:
            raw = image

        # create a dictionary to store some of the intermediate images
        intermediates = {}
        intermediates['original'] = raw


        # remove cosmic rays (now moved to an earlier stage)
        image = raw

        # define an image marking the distance from the ridge of the trace
        distanceFromTrace = np.abs(self.s - self.traceCenter(self.w))

        # define the extraction aperture (Gaussian profile for wavelength extraction, boxcar for stellar flux)
        self.speak('defining the extraction mask')
        if arc:
          intermediates['extractMask'] = self.images['RoughLSF']
          subtractsky = False
          self.speak('using a Gaussian approximation to the line-spread function (for arc extraction)', 2)
        else:
          intermediates['extractMask'] = self.ones()
          self.speak('using a boxcar', 2)

        # replaced self.obs.nFWHM*width with "obs.extractionWidth"
        intermediates['extractMask'][distanceFromTrace > self.obs.extractionWidth] = 0

        # load the (custom-defined) sky estimation mask
        intermediates['skyMask'] = self.images['skyMask']

        # wavelength calibrate the spectrum, if you can
        try:
          self.extracted['wavelength'] = self.wavelengthcalibrator.pixelstowavelengths(self.waxis)
        except AttributeError:
          self.extracted['wavelength'] = None


        # keep track of the cosmics that were rejected along the important columns
        try:
            intermediates['smearedcosmics'] = self.mask.ccd.cosmicdiagnostic[self.xstart:self.xend].reshape(1,self.xend - self.xstart)*np.ones_like(image)/(self.yend - self.ystart).astype(np.float)
            self.extracted['cosmicdiagnostic'] = (intermediates['smearedcosmics']*intermediates['extractMask']).sum(self.sindex)
            self.speak('the cosmic over-correction diagnostic is {0}'.format(np.sum(self.extracted['cosmicdiagnostics'])))
        except:
            self.speak("couldn't find any cosmic over-correction diagnostics for this frame")

        # subtract the sky, if requested
        if subtractsky:
            # estimate a 1D sky spectrum by summing over the masked region
            self.speak('estimating a sky background image')

            # currently taking median of a masked array -- would be better to fit low-order polynomial to each column
            ma = np.ma.MaskedArray(image*intermediates['skyMask']/self.images['NormalizedFlat'], intermediates['skyMask']==0)
            skyperpixel = np.ma.median(ma, self.sindex).data

            # store the sky estimation (in both 1D and 2D)
            self.extracted['sky_median'] = skyperpixel*intermediates['extractMask'].sum(self.sindex)
            intermediates['sky_median'] = np.ones_like(image)*skyperpixel.reshape((self.waxis.shape[0],1))

            # try a better sky estimation
            #self.display.one(image*intermediates['skyMask']/self.images['NormalizedFlat'])
            intermediates['sky'] = zachopy.twod.polyInterpolate(image/self.images['NormalizedFlat'], intermediates['skyMask'] == 0, order=2, visualize=False)
            self.extracted['sky'] = (intermediates['sky']*intermediates['extractMask']).sum(self.sindex)

            # for raw counts, weight by extraction mask, divide by flat, subtract the sky
            self.extracted['raw_counts'] = (intermediates['extractMask']*image/self.images['NormalizedFlat']).sum(self.sindex) - self.extracted['sky']

            # for testing, save a non-flatfielded version extraction, just to make sure
            self.extracted['no_flat'] =  (intermediates['extractMask']*(image - intermediates['sky']*self.images['NormalizedFlat'])).sum(self.sindex)


            # store the 2D sky subtracted image
            intermediates['subtracted'] = image/self.images['NormalizedFlat'] - intermediates['sky']

            if self.obs.slow:
                writeFitsData(intermediates['subtracted'], self.extractedFilename.replace('extracted', 'subtracted').replace('npy', 'fits'))

            # store a few more diagnostics
            self.extracted['centroid'] = np.nansum(self.s*intermediates['subtracted']*intermediates['extractMask'],self.sindex)/np.nansum(intermediates['subtracted']*intermediates['extractMask'], self.sindex)
            self.extracted['width'] = np.sqrt(np.nansum(self.s**2*intermediates['subtracted']*intermediates['extractMask'], self.sindex)/np.nansum(intermediates['subtracted']*intermediates['extractMask'], self.sindex) - self.extracted['centroid']**2)
            self.extracted['peak'] = np.max(intermediates['extractMask']*image, self.sindex)
            #if self.visualize:
            #    self.plot()

        else:
            self.extracted['raw_counts'] = np.nansum(intermediates['extractMask']*image/self.images['NormalizedFlat'], self.sindex)


        '''if True:
            self.figure = plt.figure(figsize=(1.0, ))
                                self.aximage.imshow(np.transpose(np.log(self.images['Science'])), cmap='gray', \
                                              extent=extent, \
                                              interpolation='nearest', aspect='auto', \
                                              vmin=np.log(1), vmax=np.log(np.nanmax(self.images['Science'])*0.01))
                                self.aximage.imshow(np.transpose(mask), alpha=0.1, cmap='winter_r', \
                                            extent=extent, \
                                            interpolation='nearest', aspect='auto')
        '''
        if False:#(arc == False):
            shape = intermediates['original'].shape
            aspect = shape[0]/(self.obs.extractionWidth*6.0)

            plt.ioff()
            dpi = 100
            plt.figure('extraction', figsize=(shape[0]*1.0/dpi,self.obs.extractionWidth*6.0/dpi), dpi=dpi)
            gs = plt.matplotlib.gridspec.GridSpec(2,1,left=0, right=1, bottom=0, top=1, wspace=0, hspace=0)
            axorig = plt.subplot(gs[0])
            axsubt = plt.subplot(gs[1], sharex=axorig, sharey=axorig)
            plt.cla()

            def scale(x):
                l = np.log(x)
                l[np.isfinite(l) == False] = 0.0
                return l

            axorig.imshow(np.transpose(scale(intermediates['original'])), interpolation='nearest', cmap='gray', aspect='equal', vmin=np.log(1), vmax=np.log(50000.0), extent=[self.waxis.min(), self.waxis.max(), self.saxis.min(), self.saxis.max()])


            axsubt.imshow(np.transpose(scale(intermediates['subtracted'])), interpolation='nearest', cmap='gray', aspect='equal', vmin=np.log(1), vmax=np.log(50000.0), extent=[self.waxis.min(), self.waxis.max(), self.saxis.min(), self.saxis.max()])

            axorig.set_xlim(self.waxis.min(), self.waxis.max())
            axorig.set_ylim((self.traceCenter(self.waxis) - 1.5*self.obs.extractionWidth).min(), (self.traceCenter(self.waxis) + 1.5*self.obs.extractionWidth).max())

            for ax in [axorig, axsubt]:
                plt.sca(ax)
                plt.plot(self.waxis, self.traceCenter(self.waxis), alpha=0.5, linewidth=3, color='seagreen')
                plt.plot(self.waxis, self.traceCenter(self.waxis) + self.obs.extractionWidth, alpha=1, linewidth=3, color='springgreen')
                plt.plot(self.waxis, self.traceCenter(self.waxis) - self.obs.extractionWidth, alpha=1, linewidth=3, color='springgreen')
                plt.setp(ax.get_xticklabels(), visible=False)
                plt.setp(ax.get_yticklabels(), visible=False)

            figfilename = self.extractedFilename.replace('extracted', 'image').replace('npy', 'png')
            plt.savefig(figfilename)
            self.speak('saved image to {0}'.format(figfilename))


            # display the calibration stamps
            '''
            self.display.window.set("frame 0")
            self.display.rgb(intermediates['original'], intermediates['extractMask'], intermediates['skyMask'])
            self.displayStamps(intermediates)
            answer = 'y'#self.input('Do you like the extraction intermediates for aperture {0}?'.format( self.name)).lower()
            assert("n" not in answer)
            if 'y' in answer:
                self.visualize=False
            '''

        # KLUDGE to prevent turnover in the wavelength calibration
        #w = self.extracted['wavelength']
        #slopeofwavelength = w[1:] - [:-1]

        #ok = slope

        #np.save(self.extractedFilename.replace('self.extracted', 'intermediates'), intermediates)
        np.save(self.extractedFilename, self.extracted)
        self.speak('saved extracted spectrum to {0}'.format(self.extractedFilename))

        if arc == False:
            #    self.visualize=True
            self.interpolate(remake=remake)

    return self.extracted

  def recalibrate(self, n):
      '''just redo the wavelength calibration'''
      self.n = n
      self.speak('recalibrating wavelengths for {0}'.format(self.n))
      # load the spectrum
      self.extracted = np.load(self.extractedFilename)[()]

      # recalibrate the spectrum
      self.extracted['wavelength'] = self.wavelengthcalibrator.pixelstowavelengths(self.waxis)

      # resave the new spectrum
      np.save(self.extractedFilename, self.extracted)

      # redo the interpolation
      self.interpolate(remake=True)

  def interpolate(self, remake=False):
        '''Interpolate the spectra onto a common (uniform) wavelength scale.'''


        # decide which values to supersample
        self.keys = ['sky',  'centroid', 'width', 'peak', 'raw_counts']
        self.limits = {}

        #
        self.supersampledFilename = self.extractedFilename.replace('extracted', 'supersampled')

        try:
            assert(remake == False)
            self.supersampled = np.load(self.supersampledFilename)
            self.speak('loaded supersampled spectrum from {0}'.format(self.extractedFilename))
        except:

            # define an empty cube that we're going to populate with spectra
            if self.visualize:
                # set up a plot window to show how the interpolation is going
                plt.figure('interpolating spectra')
                self.ax_supersampled = {}
                sharex=None
                gs = plt.matplotlib.gridspec.GridSpec(len(self.keys),1,hspace=0,wspace=0)

            # pull out the extracted wavelength
            wavelength = self.extracted['wavelength']

            # set up a fine, common wavelength grid onto which everything will be interpolated
            try:
                self.supersampled['wavelength']
            except:
                if self.obs.grism == 'vph-all':
                    commonwavelength = np.arange(4000, 10500)
                elif self.obs.grism == 'vph-red':
                    commonwavelength = np.arange(5000, 10500)

                # calculate the number of pixels that go into each wavelength bin
                dw_original = wavelength[1:] - wavelength[0:-1]
                dw_new = np.ones_like(commonwavelength)
                interpolation = scipy.interpolate.interp1d(wavelength[0:-1], dw_original, bounds_error=False)
                dnpixelsdw = dw_new/interpolation(commonwavelength)
                self.supersampled = {}

                self.supersampled['wavelength'] = commonwavelength
                self.supersampled['dnpixelsdw'] = dnpixelsdw

            # loop over the measurements
            sharex=None
            for i in range(len(self.keys)):
                key = self.keys[i]

                # set up the plots
                if self.visualize:
                    self.ax_supersampled[key] = plt.subplot(gs[i], sharex=sharex)
                    sharex = self.ax_supersampled[key]

                else:
                    # clear the plots
                    if self.visualize:
                        self.ax_supersampled[key].cla()

                # supersample onto the grid
                self.supersampled[key] = zachopy.oned.supersample(wavelength, self.extracted[key], self.supersampled['wavelength'])

                # plot demonstration
                if self.visualize and False:
                    self.ax_supersampled[key].set_ylabel(key)
                    self.ax_supersampled[key].plot(wavelength, self.extracted[key], color='black', alpha=0.5)
                    self.ax_supersampled[key].plot(self.supersampled['wavelength'], self.supersampled[key], color='red', alpha=0.5)
                    self.ax_supersampled[key].set_xlim(np.min(self.supersampled['wavelength'])-200, np.max(self.supersampled['wavelength']) + 200)
                    try:
                        self.limits[key]
                    except:
                        lims = np.min(self.supersampled[key]), np.max(self.supersampled[key])
                        span = lims[1] - lims[0]
                        nudge = .2
                        self.limits[key] = np.maximum(lims[0] - span*nudge, 0),  (lims[1] + span*nudge)
                    self.ax_supersampled[key].set_ylim(*self.limits[key])

                    if key != self.keys[-1]:
                        plt.setp(self.ax_supersampled[key].get_xticklabels(), visible=False)

            if self.visualize:
                plt.draw()

            #self.input('do you like the interpolation for {0}?'.format(self.name))
            np.save(self.supersampledFilename, self.supersampled)
            self.speak('saved supersampled spectrum to {0}'.format(self.supersampledFilename))
        return self.supersampled


  def plot(self, coordinate='wavelength', sharex=True, filename=None):
    '''Plot extracted spectrum.'''

    # switch to wavelength coordinates, if necessary (and possible)
    x = self.extracted['w']
    if coordinate=='wavelength':
      x = self.extracted['wavelength']

    # select the large figure, and clear it
    keys = self.extracted.keys()
    try:
        self.eax
    except:
        self.efig = plt.figure('extracted spectrum')
        gs = plt.matplotlib.gridspec.GridSpec(len(keys), 1, wspace=0, hspace=0)
        self.eax = []
        sharex = None
        for i in range(len(keys)):
            self.eax.append(plt.subplot(gs[i], sharex=sharex))
            sharex = self.eax[0]

    for i in range(len(keys)):
      self.eax[i].cla()
      self.eax[i].plot(x,self.extracted[keys[i]])
      self.eax[i].set_ylabel(keys[i])
    self.eax[-1].set_xlim(x.min(), x.max())
    self.eax[-1].set_xlabel(coordinate)
    self.eax[0].set_title(self.name)
    plt.draw()
    if filename != None:
      thisfig = plt.gcf()
      thisfig.savefig(filename, format='png')

  def removeCosmics(self, stamp):
    '''Subtract out cosmic rays from the science stamp.'''
    # correct for cosmic rays
    outliers = (stamp - self.images['Science'])/self.images['ScienceStdDev'] > 5
    # a bit of a kludge!
    stamp[outliers] = self.images['Science'][outliers]
    return stamp
