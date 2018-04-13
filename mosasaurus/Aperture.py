from .imports import *
from .Tools import *
from .Trace import Trace
from .WavelengthCalibrator import WavelengthCalibrator
from craftroom.cmaps import one2another
from craftroom.displays.movie import Movie
ignorekw = dict(cmap=one2another('palevioletred', 'palevioletred', alphatop=0.75, alphabottom=0.0),
                    interpolation='nearest',
                    aspect='auto',
                    zorder=1,
                    origin='lower',
                    vmin=0.0, vmax=1.0,)

extractedgekw = dict(linewidth=4, color='darkorange', alpha=0.75, zorder=2)
extractmaskkw = dict(cmap=one2another('aquamarine', 'aquamarine', alphabottom=0.0, alphatop=1.0),
                        interpolation='nearest',
                        aspect='auto',
                        zorder=100,
                        vmin=0.0, vmax=1.0,
                        origin='lower')
def imscale(x):
    l = np.log(x)
    l[np.isfinite(l) == False] = np.min(l[np.isfinite(l)]) - 1.0
    return l


plt.ion()
class Aperture(Talker):
  '''Aperture objects handle individual apertures (e.g. slits),
        storing data and handling extraction. '''

  def __init__(self, x, y, mask, **kwargs):
    '''Initialize an aperture, taking a position on the detector as input.'''

    # decide whether or not this Reducer is chatty
    Talker.__init__(self, **kwargs)

    self.mask = mask
    self.calib = self.mask.calib
    self.obs = self.calib.obs
    self.instrument = self.obs.instrument
    self.display = self.calib.display
    #zachods9('aperture', xsize=800, ysize=200, rotate=90)
    self.setup(x,y)
    self.createCalibStamps()
    # testing!
    #self.trace.fit()
    #self.createTrace()
    #self.createWavelengthCal()

  #@property
  #def directory(self):
      # something of a KLUDGE
      #return self.mask.reducer.extractionDirectory

  @property
  def visualize(self):
    '''should we visualize the steps in this process?'''
    return self.mask.reducer.visualize

  def setup(self,x,y):
    '''Setup the basic geometry of the aperture.'''

    if self.instrument.name == 'LDSS3C':
      self.x = x
      self.y = y
      self.maskWidth = self.instrument.extractiondefaults['spatialsubarray']
      #(self.obs.skyWidth + self.obs.skyGap)*2 #+20 +
      blueward = self.instrument.extractiondefaults['stampwavelengthblueward']
      redward = self.instrument.extractiondefaults['stampwavelengthredward']
      self.ystart = np.maximum(y - blueward, 0).astype(np.int)
      self.yend = np.minimum(y + redward, self.instrument.ysize).astype(np.int)
      self.xstart = np.maximum(x - self.maskWidth, 0).astype(np.int)
      self.xend = np.minimum(x + self.maskWidth, self.instrument.xsize).astype(np.int)
      # remember python indexes arrays by [row,column], which is opposite [x,y]
      x_fullframe, y_fullframe = np.meshgrid(np.arange(self.calib.images['science'].shape[1]),
                        np.arange(self.calib.images['science'].shape[0]))
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
      self.directory = os.path.join(self.mask.reducer.extractionDirectory, self.name)
      mkdir(self.directory)
      self.speak("created a spectroscopic aperture at ({0:.1f}, {1:.1f})".format(self.x, self.y))

    elif self.obs.instrument == 'IMACS':
      self.x = x
      self.y = y
      self.maskWidth = self.obs.subarray
      #(self.obs.skyWidth + self.obs.skyGap)*2 #+20 +
      self.ystart = np.maximum(y - self.obs.blueward, 0).astype(np.int)
      self.yend = np.minimum(y + self.obs.redward, self.obs.ysize).astype(np.int)
      self.xstart = np.maximum(x - self.maskWidth, 0).astype(np.int)
      self.xend = np.minimum(x + self.maskWidth, self.obs.xsize).astype(np.int)
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
    calibfilename = os.path.join(self.directory, 'calibStamps_{0}.npy'.format(self.name))
    try:
      self.images = np.load(calibfilename)[()] # i don't understand why i need the "empty tuple" but the internet says so
      self.speak("loaded calibration stamps from {0}".format(calibfilename))
    except:
      self.speak("cutting them for the first time out of the stitched master images")

      # define an empty dictionary to store the calibration stamps
      self.images = {}

      # include the coordinate system over the grid of the image
      self.images['s'] = self.s
      self.images['w'] = self.w
      self.speak("populating calibration stamps with the (s)patial and (w)avelength pixel coordinate images")


      # cut out stamps from the big images
      interesting = ['science' ,
                    'flat',
                    'He', 'Ne', 'Ar',
                    'badpixels', 'dark', 'bias']
      for k in interesting:
        self.images[k] = self.stamp(self.calib.images[k])

      #self.speak('these are the stamps before interpolating over bad pixels')
      #self.displayStamps(self.images, keys = ['science', 'flat', 'badpixels'])
      #self.input('', prompt='(press return to continue)')
      for k in interesting:
          if k != 'badpixels':
              self.images[k] = craftroom.twod.interpolateOverBadPixels(self.images[k], self.images['badpixels'])

      #self.speak('and these are they after interpolating over bad pixels')
      #self.displayStamps(self.images, keys = ['science', 'flat', 'badpixels'])
      #self.input('', prompt='(press return to continue)')

      # subtract dark from everything but the dark
      #for k in self.images.keys():
      #	if k is not 'dark':
      #		self.images[k] -= self.images['dark']

      # this is a rough flat - will be refined later in create NormalizedFlat
      self.images['RoughFlat'] = self.images['flat']/np.median(self.images['flat'], self.sindex).reshape(self.waxis.shape[0], 1)

      np.save(calibfilename, self.images)
      self.speak("saved calibration stamps to {0}".format(calibfilename))

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
    tracefilename = os.path.join(self.directory, 'trace_{0}.npy'.format(self.name))
    skyfilename = os.path.join(self.directory, 'skyMask_{0}.npy'.format(self.name))

    # define the trace object (will either load saved, or remake)
    self.trace = Trace(self)


    self.images['RoughLSF'] = np.exp(-0.5*((self.s - self.trace.traceCenter(self.w))/self.trace.tracefitwidth)**2)
    #self.displayTrace()

    # se if normalized flat exists in self.images, if not make one!
    try:
      self.images['NormalizedFlat']
    except(KeyError):
      self.speak("creating normalized flat for {0}".format(self.name))
      self.createNormalizedFlat()

  def displayTrace(self):
      for i, width in enumerate(self.trace.extractionwidths):
          self.display.new(frame=i)
          self.display.rgb( self.images['science'],
                            self.trace.extractionmask(width),
                            self.trace.skymask(width))

  def createNormalizedFlat(self, visualize=True):

      # create a normalized flat field stamp, dividing out the blaze + spectrum of quartz lamp
      # make sure the normalization only happens over regions of interest - don't want to include parts with no exposure

      normfilename = os.path.join(self.directory, 'normFlat_{0}.pdf'.format(self.name))
      calibfilename = os.path.join(self.directory, 'calibStamps_{0}.npy'.format(self.name))

      #testing using no flat
      #self.images['NormalizedFlat'] = np.ones(self.images['flat'].shape)

      # use the skymask function to get a sky mask, not aperture-width based
      infmask = self.trace.skymask(-np.inf)
      normmask = np.zeros_like(self.images['science'])

      # figure out the sky boundaries so everything inbetween can be set to True and a median can be taken across the right indies in each row
      minarray, maxarray = np.zeros(self.waxis.shape[0]), np.zeros(self.waxis.shape[0])
      medfilt_1d = np.zeros(self.waxis.shape[0])
      for i in range(self.waxis.shape[0]):
          whereone = np.where(infmask[i])
          minindex, maxindex = np.min(whereone), np.max(whereone)
          minarray[i], maxarray[i] = minindex, maxindex
          normmask[i][minindex:maxindex+1] = 1
          medfilt_1d[i] = np.median(self.images['flat'][i][minindex:maxindex+1])

      # divide masked flat by the 1d_median filter, reshaped for division purposes
      self.images['NormalizedFlat'] = (self.images['flat']*normmask)/medfilt_1d.reshape(self.waxis.shape[0], 1) 
      # need values outside of extraction to be 1, not 0, for future division purposes
      self.images['NormalizedFlat'] += np.abs(normmask - 1)

      # visualize and save NormalizedFlat
      plt.figure('normalized flat')
      ax = plt.subplot()
      ax.imshow(self.images['NormalizedFlat'].T, cmap='gray', aspect='auto', vmin=0.95, vmax=1.05)
      ax.set_xlabel('pixels')
      ax.set_ylabel('pixels')
      if visualize:
          # flat field has been medianed so most values should center on 1; want to be able to see finer-scale structure
          self.display.one(self.images['NormalizedFlat'], aspect='auto', vmin=0.95, vmax=1.05) 
          self.display.run()
          answer = self.input("Did you like the NormalizedFlat for this stamp? [Y,n]").lower()
      assert('n' not in answer)

      plt.savefig(normfilename)
      self.speak("saved normalized flat to {0}".format(normfilename))

      np.save(calibfilename, self.images)
      self.speak("saved normalized flat to calibration stamp {0}".format(calibfilename))

  def createWavelengthCal(self, remake=False):
    '''Populate the wavelength calibration for this aperture.'''
    self.speak("populating wavelength calibration")
    self.wavelengthcalibrator = WavelengthCalibrator(self)

  def ones(self):
    '''Create a blank array of ones to fill this aperture.'''
    return np.ones_like(self.images['science'])

  @property
  def extractedFilename(self):
      return os.path.join(self.directory, 'extracted_{}.npy'.format(self.exposureprefix))

  def loadExtracted(self, n):
    self.exposureprefix = n
    self.extracted = np.load(self.extractedFilename)[()]
    self.speak('loaded extracted spectrum from {0}'.format(self.extractedFilename))

  def extract(self, n=0, image=None, subtractsky=True, arc=False, remake=False):
    '''Extract the spectrum from this aperture.'''

    # make sure this aperture is pointed at the desired image
    self.exposureprefix = n
    try:
        # try to load a previously saved spectrum
        assert(remake == False)
        assert(arc == False)
        assert(os.path.exists(self.extractedFilename))

        #self.extracted = np.load(self.extractedFilename)[()]
        #self.speak('raw extracted file {0} already exists'.format(self.supersampledFilename))
        return
    except (AssertionError, IOError):
        self.speak('extracting spectrum from image {}'.format(self.exposureprefix))
        # make sure the trace is defined
        try:
            self.trace.traceCenter
        except AttributeError:
            self.createTrace()


        # reassign the image number for this, because createwavelength cal probably reset it
        self.exposureprefix = n

        # make sure we don't subtract the sky, if the arcs are set
        subtractsky = subtractsky & (arc == False)


        # create a dictionary to store the self.extracted spectra
        self.extracted = {}
        self.extracted['w'] = self.waxis

        if image is None:
            # make sure the right data have been loaded
            #assert(self.mask.ccd.exposureprefix == n)
            self.mask.load(n)

            # extract the raw science stamp from this mask image
            raw = self.stamp(self.mask.ccd.data)
        else:
            raw = image

        # create a dictionary to store some of the intermediate images
        self.intermediates = {}
        self.intermediates['original'] = raw

        # remove cosmic rays (now moved to an earlier stage)
        image = raw



        # define the extraction aperture (Gaussian profile for wavelength extraction, boxcar for stellar flux)
        for i, width in enumerate(self.trace.extractionwidths):
            # create a dictionary to store this extraction width
            self.extracted[width], self.intermediates[width] = {}, {}
            self.speak('extracting spectrum for width of {}'.format(width))
            self.intermediates[width]['extractMask'] = self.trace.extractionmask(width)
            if arc:
                self.intermediates[width]['extractMask'] = self.intermediates[width]['extractMask']*self.images['RoughLSF']
                self.speak('using a Gaussian approximation to the line-spread function (for arc extraction)', 2)

            # load the appropriate sky estimation mask
            self.intermediates[width]['skyMask'] = self.trace.skymask(width)

            # keep track of the cosmics that were rejected along the important columns
            if self.instrument.zapcosmics:
                try:
                    if arc == False:
                        self.intermediates[width]['smearedcosmics'] = self.mask.ccd.cosmicdiagnostic[self.xstart:self.xend].reshape(1,np.round(self.xend - self.xstart).astype(np.int))*np.ones_like(image)/(self.yend - self.ystart).astype(np.float)
                        self.extracted[width]['cosmicdiagnostic'] = (self.intermediates[width]['smearedcosmics']*self.intermediates[width]['extractMask']).sum(self.sindex)
                        self.speak('the cosmic over-correction diagnostic is {0}'.format(np.sum(self.extracted[width]['cosmicdiagnostic'])))
                except (IOError, AttributeError):
                    self.speak("couldn't find any cosmic over-correction diagnostics for this frame")

            #import sys
            #sys.exit("Breaking here. Check it out.")


            # subtract the sky, if requested
            if subtractsky:

                # estimate a 1D sky spectrum by summing over the masked region
                self.speak('estimating a sky background image')
                #mask = np.ones(self.images['BadPixels'].shape)

                '''
                #####################################################
                # hzdl modification: making extraction and sky windows based on FWHM of each row in the cross-dispersion direction
                mask = (self.images['BadPixels']-1)*-1
                subimage = (image/self.images['NormalizedFlat']*mask)#[:,boxcuts[0]:boxcuts[-1]]
                disp, crossdisp = subimage.shape
                FWHM = np.zeros(disp)
                edge1, edge2 = np.zeros(disp), np.zeros(disp)

                for i in range(disp):
                    
                    # make a psf by cross-dispersion column; multiply by extraction mask so neighboring spectra are not included
                    psf = subimage[i]*self.intermediates[width]['extractMask'][i]
                    # create a spline to fit to the psf; this will give the FWHM roots
                    spline = spi.UnivariateSpline(range(len(psf)), psf-np.max(psf)/2.)
                    r1, r2 = spline.roots() # find the roots

                    FWHM[i] = r2-r1
                    edge1[i] = r1
                    edge2[i] = r2

                # up to here edges will be the same for every extraction width;
                # extend edges out by width/2 times FWHM/2 on each side
                edge1 -= (FWHM/2.)*(width/2.)
                edge2 += (FWHM/2.)*(width/2.)

                # smooth edges
                edge1spline = spi.UnivariateSpline(range(len(edge1)), edge1)
                edge1smooth = edge1spline(range(len(edge1)))
                edge2spline = spi.UnivariateSpline(range(len(edge2)), edge2)
                edge2smooth = edge2spline(range(len(edge2)))

                # this will become the new extraction mask; make sure to include partial pixels 
                newextractmask = np.zeros(subimage.shape)
                for i in range(disp):
                    # hard edges will have values of 1; soft edges have partial pixel values; be careful of indexing
                    hardedge1, hardedge2 = int(np.ceil(edge1smooth[i])), int(np.floor(edge2smooth[i]))
                    softedge1, softedge2 = int(np.floor(edge1smooth[i])), int(np.floor(edge2smooth[i]))
                    newextractmask[i][hardedge1:hardedge2] = 1
                    newextractmask[i][softedge1] = hardedge1 - edge1smooth[i]
                    newextractmask[i][softedge2] = edge2smooth[i] - hardedge2

                # change sky mask so that there are no overlaps with new extract mask
                newskymask = self.intermediates[width]['skyMask']
                # find where the new extraction mask is overlapping with the sky mask
                removesky = newextractmask + newskymask > 1
                # remove the overlapped values from the sky mask
                newskymask -= removesky

                # replace sky mask and extraction mask with new versions
                self.intermediates[width]['extractMask'] = newextractmask
                self.intermediates[width]['skyMask'] = newskymask

                # end hzdl modification
                ###############################
                '''

                # do polynomial interpolation along each column to estimate sky
                self.intermediates[width]['sky'] = craftroom.twod.polyInterpolate(image/self.images['NormalizedFlat'], self.intermediates[width]['skyMask'] == 0, order=2, visualize=False)

                self.extracted[width]['sky'] = (self.intermediates[width]['sky']*self.intermediates[width]['extractMask']).sum(self.sindex)

                # for raw counts, weight by extraction mask, divide by flat, subtract the sky
                self.extracted[width]['raw_counts'] = (self.intermediates[width]['extractMask']*image/self.images['NormalizedFlat']).sum(self.sindex) - self.extracted[width]['sky']

                # for testing, save a non-flatfielded version extraction, just to make sure
                self.extracted[width]['no_flat'] =  (self.intermediates[width]['extractMask']*(image - self.intermediates[width]['sky']*self.images['NormalizedFlat'])).sum(self.sindex)

                # store the 2D sky subtracted image
                self.intermediates[width]['subtracted'] = image/self.images['NormalizedFlat'] - self.intermediates[width]['sky']

                '''
                # optimal extraction
                # define the parameters that will act as inputs to the optimal extraction code optspex.py
                subdata = image/self.images['NormalizedFlat'] - self.intermediates[width]['sky']
                # draw a box around the spectral trace; box must be big enough to include the whole trace, but not so big as to include areas that were not sky-subtracted
                boxcuts = []
                for i in range(subdata.shape[1]):
                    if 1. in self.intermediates[width]['extractMask'][:,i]: boxcuts.append(i)
                subdata = subdata[:, boxcuts[0]:boxcuts[-1]]
                submask = (((self.images['BadPixels']-1)*-1)*mask)[:, boxcuts[0]:boxcuts[-1]]
                bg = intermediate_sky[:, boxcuts[0]:boxcuts[-1]]
                spectrum = self.extracted[width]['raw_counts']

                #startrace = trace.calc_trace(subdata, submask)#, nknots=750, deg=4)
                #variance = abs(subdata)
                #spec_width = 4
                #trdata, trmask, trvar, trbg = trace.realign(subdata, submask, variance, bg, spec_width, startrace)
                #trspec, trspecunc, trnewmask = OptimalExtraction.optimize(trdata.T, trmask.T, trbg.T, spectrum, 1, 0, p5thresh=5, p7thresh=5, fittype='smooth', window_len=11)

                spec, specunc, newmask = OptimalExtraction.optimize(subdata.T, submask.T, bg.T, spectrum, 1, 0, p5thresh=10, p7thresh=10, fittype='smooth', window_len=11)
                self.extracted[width]['raw_counts_optext'] = spec
                '''
                self.extracted[width]['raw_counts_optext'] = self.extracted[width]['raw_counts']


                #if self.obs.slow:
                #    writeFitsData(self.intermediates['subtracted'], self.extractedFilename.replace('extracted', 'subtracted').replace('npy', 'fits'))

                # store a few more diagnostics
                fine = np.isfinite(self.intermediates[width]['subtracted'])

                weights = fine*self.intermediates[width]['subtracted']*self.intermediates[width]['extractMask']
                self.extracted[width]['centroid'] = np.sum(weights*self.s, axis=self.sindex)/np.sum(weights, axis=self.sindex)
                #self.extracted[width]['centroid'] = np.average(
                #                        self.s,
                #                        axis=self.sindex,
                #                        weights=fine*self.intermediates[width]['subtracted']*self.intermediates[width]['extractMask'])

                reshapedFWC = self.extracted[width]['centroid'][:,np.newaxis]
                weights = fine*self.intermediates[width]['subtracted']*self.intermediates[width]['extractMask']
                weights = np.maximum(weights, 0)
                weights += 1.0/np.sum(weights)
                self.extracted[width]['width'] = np.sqrt(
                    np.average((self.s-reshapedFWC)**2, axis=self.sindex, weights=weights))

                self.extracted[width]['peak'] = np.max(self.intermediates[width]['extractMask']*image, self.sindex)

                # KLUDGE -- am I definitely catching this error later?
                #for k in ['centroid', 'width', 'peak']:
                #    assert(np.isfinite(self.extracted[width][k]).all())
            else:
                self.extracted[width]['raw_counts'] = np.nansum(self.intermediates[width]['extractMask']*image/self.images['NormalizedFlat'], self.sindex)
                # this is a kludge, to make the plotting look better for arcs
                self.intermediates[width]['sky']  = np.zeros_like(self.intermediates['original'])# + np.percentile(self.extracted[width]['raw_counts'] , 1)

        if arc==False:
            self.visualizeExtraction(width)


        # save the extracted spectra for all apertures
        np.save(self.extractedFilename, self.extracted)
        self.speak('saved extracted spectra to {0}'.format(self.extractedFilename))

    # create a supersampled version of the extracted spectra
    #if False:#arc == False:
    #    self.interpolate(remake=remake)

    # return the extracted spectrum
    return self.extracted

  def setupVisualization(self):
        self.thingstoplot = ['raw_counts', 'raw_counts_optext']#['sky', 'width',  'raw_counts']
        height_ratios = np.ones(len(self.thingstoplot) + 2)
        suptitletext = '{}, {}'.format(self.name, self.exposureprefix)
        try:
            self.suptitle.set_text(suptitletext)
            self.plotted
        except AttributeError:
            self.figure = plt.figure(figsize=(20,12), dpi=50)
            gs = plt.matplotlib.gridspec.GridSpec(len(self.thingstoplot) + 2, self.trace.numberofapertures,
                                                    hspace=0.1, wspace=0.05,
                                                    height_ratios = height_ratios )
            self.suptitle = self.figure.suptitle(suptitletext, fontsize=17)
            # create all the axes
            self.ax = {}
            self.plotted = {}
            sharex = None
            shareapertures, sharesubtracted = None, None
            sharey = {}
            for thing in self.thingstoplot:
                sharey[thing] = None
            for i, width in enumerate(self.trace.extractionwidths):
                self.ax[width], self.plotted[width] = {}, {}
                self.ax[width]['apertures'] = plt.subplot(gs[-2, i], sharex=sharex, sharey=shareapertures)
                sharex = self.ax[width]['apertures']
                shareapertures=self.ax[width]['apertures']
                self.ax[width]['subtracted'] = plt.subplot(gs[-1, i], sharex=sharex, sharey=sharesubtracted)
                sharesubtracted = self.ax[width]['subtracted']
                self.ax[width]['subtracted'].set_xlabel('Pixels')
                self.ax[width]['apertures'].set_ylim(*craftroom.oned.minmax(self.saxis))
                if i == 0:
                    self.ax[width]['apertures'].set_ylabel('extraction apertures')
                    self.ax[width]['subtracted'].set_ylabel('sky-subtracted, and coarsely rectified')
                else:
                    plt.setp(self.ax[width]['subtracted'].get_yticklabels(), visible=False)
                    plt.setp(self.ax[width]['apertures'].get_yticklabels(), visible=False)
                plt.setp(self.ax[width]['apertures'].get_xticklabels(), visible=False)

                for j, thing in enumerate(self.thingstoplot):
                    self.ax[width][thing] = plt.subplot(gs[j,i], sharex=sharex, sharey=sharey[thing])
                    self.ax[width][thing].locator_params(nbins=3)
                    sharey[thing] = self.ax[width][thing]

                    if j == 0:
                        self.ax[width][thing].set_title('{} pix.'.format(width))
                    if i == 0:
                        self.ax[width][thing].set_ylabel(thing)
                    else:
                        plt.setp(self.ax[width][thing].get_yticklabels(), visible=False)
                    plt.setp(self.ax[width][thing].get_xticklabels(), visible=False)

                    if thing == 'width':
                        self.ax[width][thing].set_ylim(0,np.max(self.trace.extractionwidths)/3.5)
                    if thing == 'centroid':
                        self.ax[width][thing].set_ylim(np.min(self.trace.traceCenter(self.waxis))-5, np.max(self.trace.traceCenter(self.waxis))+5)
                    if thing == 'raw_counts':
                        self.ax[width][thing].set_ylim(0, np.percentile(self.extracted[width]['raw_counts'], 99)*1.5)
                    if thing == 'raw_counts_optext':
                        self.ax[width][thing].set_ylim(0, np.percentile(self.extracted[width]['raw_counts_optext'], 99)*1.5)

  def visualizeExtraction(self, width):
        # make sure the axes have been setup
        self.setupVisualization()

        widths = [k for k in self.extracted.keys() if type(k) != str]
        for width in widths:
            try:
                self.plotRectified(width)
            except KeyError:
                self.speak('skipping sky-subtracted plot')
            self.plotApertures(width)
            self.plotExtracted(width)

        self.figure.savefig(self.extractedFilename.replace('npy', 'pdf'))

        '''
        if 'raw_counts_optext' in self.extracted[width]:
            self.setupVisualization('raw_counts_optext')

            for width in widths:
                try:
                    self.plotRectified(width)
                except KeyError:
                    self.speak('skipping sky-subtracted plot')
                self.plotApertures(width)
                self.plotExtracted(width)

            savename = self.extractedFilename[:-4]+'optext.pdf'
            self.figure.savefig(savename)
        '''
  @property
  def widths(self):
    '''list all the widths available for this extraction'''
    return np.array([k for k in self.extracted.keys() if type(k) != str])

  def plotRectified(self, width, ax=None):
        '''plot the rectified, sky-subtracted image'''
        if ax is None:
            ax = self.ax[width]['subtracted']

        # create interpolators to get indices from coordinates
        #wtoindex = scipy.interpolate.interp1d(self.waxis, np.arange(len(self.waxis)))
        #stoindex = scipy.interpolate.interp1d(self.saxis, np.arange(len(self.saxis)))

        # create grid of coordinates
        ok = self.intermediates[width]['skyMask'] > 0
        offsets = self.s[ok] - self.trace.traceCenter(self.w[ok])
        ylim = craftroom.oned.minmax(offsets)

        rectifiedw, rectifieds = np.meshgrid(self.waxis, np.arange(*ylim))
        rectifieds += self.trace.traceCenter(rectifiedw)

        indexw = np.interp(rectifiedw, self.waxis, np.arange(len(self.waxis))).astype(np.int)
        indexs = np.interp(rectifieds, self.saxis, np.arange(len(self.saxis))).astype(np.int)
        #indexw, indexs = wtoindex(rectifiedw).astype(np.int), stoindex(rectifieds).astype(np.int)
        # the s coordinate might be outside the size of the subarray
        indexs = np.maximum(np.minimum(indexs, len(self.saxis)-1), 0)
        rectified = self.intermediates[width]['subtracted'][indexw, indexs]#.reshape(rectifiedw.shape)
        ignore = ((self.intermediates[width]['skyMask'] + self.intermediates[width]['extractMask']) == 0)[indexw, indexs]
        ignore += (indexs == 0)
        ignore += (indexs == len(self.saxis)-1)
        ignore = ignore > 0

        #interpolator = scipy.interpolate.RectBivariateSpline(self.waxis, self.saxis, self.intermediates[width]['subtracted'])

        #rectifiedw, rectifieds = np.meshgrid(self.waxis, )
        #rectified = interpolator(rectifiedw, self.trace.traceCenter(rectifiedw) + rectifieds)
        #x, y = rectifiedw, self.trace.traceCenter(rectifiedw) + rectifieds
        #rectified = np.zeros_like(x)
        # why do I have to do this?
        #for i in range(y.shape[1]):
        #    rectified[:,i] = interpolator(x[:,i], y[:,i].T)



        # try updating an existing plot
        try:
            self.plotted[width]['subtractedandrectified'].set_data(rectified)
        except KeyError:
            nbackgrounds = 3
            vmin = -nbackgrounds*np.min([1.48*craftroom.oned.mad(self.intermediates[onewidth]['subtracted'][self.intermediates[width]['skyMask'] > 0]) for onewidth in self.widths])
            vmax = 0.001*np.percentile(self.extracted[width]['raw_counts'], 99)
            extent =[ rectifiedw.min(), rectifiedw.max(), ylim[0], ylim[1]]
            self.plotted[width]['subtractedandrectified'] = ax.imshow(rectified,
                     cmap='gray',
                     extent=extent,
                     interpolation='nearest',
                     aspect='auto',
                     zorder=0,
                     origin='lower',
                     vmin=vmin, vmax=vmax)

            self.plotted[width]['subtractedandrectified_toignore'] = ax.imshow(ignore,
                        extent=extent,
                        **ignorekw
                    )

            for side in [-1,1]:
                ax.axhline(side*width, **extractedgekw)



            ax.set_xlim(*craftroom.oned.minmax(rectifiedw))
            ax.set_ylim(*ylim)

            '''offsets = np.array([-1,1])*width
            x = self.waxis
            y = self.trace.traceCenter(x)[:,np.newaxis] + offsets[np.newaxis, :]
            self.plotted[thiswidth]['extractionedges'] = self.ax[thiswidth]['subtracted'].plot(x, y,
                                                linewidth=1,
                                                alpha=0.5,
                                                color='darkgreen',
                                                zorder = 1000)'''

  def plotApertures(self, width, ax=None):
        '''plot the extraction and sky apertures on top of the 2D exposure'''
        if ax is None:
            ax = self.ax[width]['apertures']

        image = self.intermediates['original']
        ignore = (self.intermediates[width]['skyMask'] + self.intermediates[width]['extractMask']) == 0

        # try updating an existing plot
        try:
            self.plotted[width]['apertures'].set_data(self.imscale(image.T))
        except KeyError:
            extent = [self.waxis.min(), self.waxis.max(), self.saxis.min(), self.saxis.max()]
            skylevel = np.percentile(self.intermediates[width]['sky'], 1)
            peaklevel = np.percentile(self.intermediates['original']*self.intermediates[width]['extractMask'], 99)
            vmin, vmax = skylevel, peaklevel, #np.percentile(imscale(image), [1,99])
            def imscale(z):
                buffer = 0.0001
                x = np.maximum((z - vmin)/(vmax - vmin), buffer)
                return np.log(x)
                #return np.sqrt(np.maximum(z - vmin, 0.0))
            self.imscale = imscale
            self.plotted[width]['apertures'] = ax.imshow(self.imscale(image.T),
                     cmap='gray',
                     extent=extent,
                     interpolation='nearest',
                     aspect='auto',
                     zorder=0,
                     origin='lower',
                     vmin=self.imscale(vmin), vmax=self.imscale(vmax))

            self.plotted[width]['apertures_toignore'] = ax.imshow(ignore.T,
                        extent=extent,
                        **ignorekw
                    )

            for side in [-1,1]:
                ax.plot(self.waxis, self.trace.traceCenter(self.waxis)[:,np.newaxis]+np.array([[-width, width]]), **extractedgekw)

  def plotExtracted(self, width, axes=None):
        '''plot the extract spectra and parameters'''

        if axes is None:
            axes = self.ax[width]

        for j, thing in enumerate(self.thingstoplot):
            if thing in self.extracted[width].keys():
                try:
                    self.plotted[width][thing][0].set_data(self.extracted['w'], self.extracted[width][thing])
                except KeyError:
                    self.plotted[width][thing] = \
                        axes[thing].plot(
                            self.extracted['w'],
                            self.extracted[width][thing],
                            linewidth=1, color='black')

  def movieExtraction(self):
      pattern = '{}/extracted*.pdf'.format(self.directory)
      directory = os.path.join(self.directory, 'animatedextraction/')
      filename= os.path.join(self.directory, 'animatedextraction.mp4')
      if os.path.exists(filename):
          self.speak('{} already exists; not remaking it.'.format(filename))
      else:
          self.speak('making a movie of the extraction, in ')
          self.movie = Movie(pattern=pattern, directory=directory, filename=filename)

  def wavelengthcalibrate(self, n):
    '''THIS STILL NEEDS TO GET FOLDED IN!'''
    # unless this is an arc, make sure there is a wavelength calibrator
    try:
      self.wavelengthcalibrator
    except AttributeError:
      if arc == False:
          self.createWavelengthCal()


    # wavelength calibrate the spectrum, if you can
    try:
      self.extracted['wavelength'] = self.wavelengthcalibrator.pixelstowavelengths(self.waxis)
    except AttributeError:
      self.extracted['wavelength'] = None



  def addWavelengthCalibration(self, exposureprefix, remake=False, shift=False):
      '''just redo the wavelength calibration'''
      # point at this CCD number
      self.exposureprefix = exposureprefix

      if shift:
          dir = os.path.join(self.directory, 'stretchedsupersampled')
          mkdir(dir)
          self.supersampledFilename = os.path.join(dir, 'stretchedsupersampled_{}.npy'.format(self.exposureprefix))
      else:
          dir = os.path.join(self.directory, 'supersampled')
          mkdir(dir)
          self.supersampledFilename = os.path.join(dir, 'supersampled_{}.npy'.format(self.exposureprefix))


      # only load and redo if supersampled doesn't exist
      if remake or not os.path.exists(self.supersampledFilename):
          self.speak('recalibrating wavelengths for {0}'.format(self.exposureprefix))
          # load the spectrum
          self.extracted = np.load(self.extractedFilename)[()]

          # addWavelengthCalibration the spectrum
          self.extracted['wavelength'] = self.wavelengthcalibrator.pixelstowavelengths(self.waxis)

          # resave the new spectrum
          np.save(self.extractedFilename, self.extracted)

          # redo the interpolation
          # self.visualize = False
          self.interpolate(shift=shift)
      else:
          self.speak('{0} already exists'.format(self.supersampledFilename))

  def interpolate(self, remake=False, shift=False):
        '''Interpolate the spectra onto a common (uniform) wavelength scale.'''


        # these are things where we care about the sum matching extracted to supersampled
        self.additivekeys = ['raw_counts', 'sky']
        # these are things where we want the individual values matching extracted to supersampled
        self.intrinsickeys = ['centroid', 'width', 'peak']
        # these are all the keys that will be supersampled
        self.keys = self.additivekeys + self.intrinsickeys

        try:
            assert(remake == False)
            self.supersampled = np.load(self.supersampledFilename)
            self.speak('loaded supersampled spectrum from {0}'.format(self.supersampledFilename))
        except (AssertionError, IOError):

            if shift:
                try:
                    self.stretches
                except AttributeError:
                    stretchfilename = os.path.join(self.mask.reducer.extractionDirectory,  'spectralstretch.npy')
                    self.stretches = np.load(stretchfilename)[()]
                    self.speak('loaded stretches from {}'.format(stretchfilename))

            # we're going to supersample multiple keys, to keep everything together
            nkeys = len(self.keys)

            # what are the extraction widths to consider?
            widths = np.sort([k for k in self.extracted.keys() if type(k) != str])
            napertures = len(widths)

            # define an empty axes that we're going to populate with spectra
            if self.visualize:
                suptitletext = '{}, {}'.format(self.name, self.exposureprefix)
                try:
                    self.suptitle_supersampled.set_text(suptitletext)
                    self.plotted_supersampled
                except AttributeError:
                    # set up a plot window to show how the interpolation is going
                    self.figure_supersampled = plt.figure('interpolating spectra', figsize=(20,10), dpi=60)

                    self.ax_supersampled = {}
                    sharex=None
                    self.suptitle_supersampled = plt.suptitle(suptitletext)

                    # make subplot for each key, and each aperture
                    gs = plt.matplotlib.gridspec.GridSpec(nkeys, napertures, hspace=0.1, wspace=0.05)
                    self.plotted_supersampled = {}

                try:
                    self.limits_supersampled
                except AttributeError:
                    # for plotting, we'll want to keep track of limits for each key
                    self.limits_supersampled = {}

            # pull out the extracted wavelength and column pixel number
            if shift:

                originalwavelength = self.extracted['wavelength']
                midpoint = self.stretches['midpoint']
                star = self.name
                coefficients = self.stretches['stretch'][star][self.exposureprefix], self.stretches['shift'][star][self.exposureprefix]
                dw = np.polyval(coefficients, originalwavelength - midpoint)
                wavelength = originalwavelength - dw
                phrase = 'dw = {:.4}x(w - {midpoint}){:+.4}'.format(*coefficients, midpoint=midpoint)
                self.speak('nudge wavelengths for {} by {}'.format(self.exposureprefix, phrase))
            else:
                wavelength = self.extracted['wavelength']
            pixelnumber = self.extracted['w']

            # set up a fine, common wavelength grid onto which everything will be interpolated
            try:
                self.supersampled['wavelength']
            except (AttributeError, KeyError):

                # pull out a reasonable common wavelength grid for this instrument
                commonwavelength = self.instrument.uniformwavelengths

                # what's the (assumed constant) dw for the uniform grid?
                scale = 1.0
                assert((np.diff(commonwavelength) == scale).all())

                # what fraction of an original pixel went into this new pixel?
                doriginaldnew = fluxconservingresample(wavelength,
                                                                        np.ones_like(pixelnumber),
                                                                        commonwavelength)


                self.supersampled = {}
                self.supersampled['wavelength'] = commonwavelength
                self.supersampled['fractionofapixel'] = doriginaldnew
            # loop over the measurement types
            sharex=None
            for i, key in enumerate(self.keys):

                # supersample onto the grid
                for j,width in enumerate(widths):

                    # make a combined key that includes both the key name, and the widthkey
                    widthkey = '{:04.1f}px'.format(width)
                    combinedkey = key + '_' + widthkey

                    # use flux-conserving resampling to put onto this grid
                    yoriginal = self.extracted[width][key]
                    isnanoriginal = np.isnan(yoriginal)

                    # this ones where the supersampled array was corrupted by nans
                    yclosetonan = fluxconservingresample(
                                        wavelength,
                                        isnanoriginal,
                                        self.supersampled['wavelength']) > 0

                    # supersample onto the uniform grid
                    ysupersampled = fluxconservingresample(
                                        wavelength,
                                        self.extracted[width][key],
                                        self.supersampled['wavelength'],
                                        treatnanas=0.0)

                    assert(np.isfinite(ysupersampled).any())

                    # turn the bad elements back to nans
                    #ysupersampled[yclosetonan] = np.nan

                    if key in self.additivekeys:
                        self.supersampled[combinedkey] = ysupersampled
                    elif key in self.intrinsickeys:
                        self.supersampled[combinedkey] = ysupersampled/self.supersampled['fractionofapixel']
                    else:
                        self.speak("Yikes! It's not clear if {} is an additive or intrinsic quantity!".format(key))

                    # set up the plots to visualize this
                    if self.visualize:
                        try:
                            # does this subplot already exist and just need to be cleared?
                            self.ax_supersampled[combinedkey]
                        except KeyError:
                            # or do we need to create it from scratch?
                            self.ax_supersampled[combinedkey] = plt.subplot(gs[i,j], sharex=sharex)
                            sharex = self.ax_supersampled[combinedkey]

                    # plot demonstration
                    if self.visualize:
                        if key in self.additivekeys:
                            ontheoriginalpixels = self.supersampled[combinedkey]/self.supersampled['fractionofapixel']
                        elif key in self.intrinsickeys:
                            ontheoriginalpixels = self.supersampled[combinedkey]
                        else:
                            ontheoriginalpixels = None

                        try:
                            # can we just updated existing plots?
                            plot_extracted, plot_supersampled = self.plotted_supersampled[combinedkey]
                            plot_extracted[0].set_data(wavelength, self.extracted[width][key])
                            plot_supersampled[0].set_data(self.supersampled['wavelength'], ontheoriginalpixels)
                        except KeyError:
                            # or do we have to make new ones?

                            # plot the original spectrum
                            if key in self.additivekeys:
                                ylabel = combinedkey + "\n(per original pixel)"
                            else:
                                ylabel = combinedkey
                            self.ax_supersampled[combinedkey].set_ylabel(ylabel)
                            plot_extracted = self.ax_supersampled[combinedkey].plot(wavelength, self.extracted[width][key], color='black', alpha=0.5)
                            plot_supersampled = self.ax_supersampled[combinedkey].plot(self.supersampled['wavelength'], ontheoriginalpixels, color='red', alpha=0.5)
                            self.plotted_supersampled[combinedkey] = [plot_extracted, plot_supersampled]
                            self.ax_supersampled[combinedkey].set_xlim(np.min(self.supersampled['wavelength'])-200, np.max(self.supersampled['wavelength']) + 200)
                            try:
                                self.limits_supersampled[combinedkey]
                            except:
                                ok = np.isfinite(ontheoriginalpixels)

                                lims = np.percentile(ontheoriginalpixels[ok], [5, 95])#np.nanmin(), np.nanmax(self.extracted[width][key])
                                span = np.abs(lims[1] - lims[0])
                                nudge = 0.5
                                self.limits_supersampled[combinedkey] = [(lims[0] - span*nudge),  (lims[1] + span*nudge)]
                                if 'centroid' not in key:
                                    self.limits_supersampled[combinedkey][0] = np.maximum(0, self.limits_supersampled[combinedkey][0])
                            self.ax_supersampled[combinedkey].set_ylim(*self.limits_supersampled[combinedkey])

                            if key == self.keys[0]:
                                plt.title(widthkey)
                            if key == self.keys[-1]:
                                plt.xlabel('Wavelength (angstroms)')
                            else:
                                plt.setp(self.ax_supersampled[combinedkey].get_xticklabels(), visible=False)
                            if width != widths[0]:
                                plt.setp(self.ax_supersampled[combinedkey].get_yticklabels(), visible=False)
                                self.ax_supersampled[combinedkey].set_ylabel('')

            if self.visualize:
                #plt.show()
                #plt.draw()
                supersamplingPlotFilename = self.supersampledFilename.replace('npy', 'pdf')
                plt.savefig(supersamplingPlotFilename)
                self.speak('saved a plot to {}'.format(supersamplingPlotFilename))
                #a = self.input('is this okay?')

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
    outliers = (stamp - self.images['science'])/self.images['ScienceStdDev'] > 5
    # a bit of a kludge!
    stamp[outliers] = self.images['science'][outliers]
    return stamp
