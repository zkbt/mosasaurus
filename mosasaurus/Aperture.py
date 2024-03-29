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

    # images *must* be row=wavelength, col=spatial
    self.windex = 0
    self.sindex = 1 - self.windex

    self.name = 'aperture_{0:.0f}_{1:.0f}'.format(self.x, self.y)
    self.directory = os.path.join(self.mask.reducer.extractionDirectory, self.name)
    mkdir(self.directory)
    self.speak("created a spectroscopic aperture at ({0:.1f}, {1:.1f})".format(self.x, self.y))

  def stamp(self, image):
    '''Return a postage stamp of an image, appropriate for this aperture.'''
    return image[self.ystart:self.yend, self.xstart:self.xend]

  def createCalibStamps(self, visualize=True, interactive=True):
    '''Populate the necessary postage stamps for the calibrations.'''
    self.speak("populating calibration stamps")
    calibfilename = os.path.join(self.directory, 'calibStamps_{0}.npy'.format(self.name))
    try:
      self.images = np.load(calibfilename, allow_pickle=True)[()] # i don't understand why i need the "empty tuple" but the internet says so
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
      #interesting = ['science', 'flat', 'He', 'Ne', 'Ar','badpixels', 'dark', 'bias']
      interesting = ['science', 'badpixels'] + self.instrument.detectorcalibrations + self.instrument.arclamps
      for k in interesting:
        self.images[k] = self.stamp(self.calib.images[k])

      #self.speak('these are the stamps before interpolating over bad pixels')
      #self.displayStamps(self.images, keys = ['science', 'flat', 'badpixels'])
      #self.input('', prompt='(press return to continue)')
      #self.images['badpixels'] = np.zeros_like(self.images['science'])
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
      try:
          self.images['RoughFlat'] = self.images['flat']/np.median(self.images['flat'], self.sindex).reshape(self.waxis.shape[0], 1)
      except(KeyError):
          self.images['flat'] = np.ones_like(self.images['science'])
          self.images['RoughFlat'] = np.ones_like(self.images['science'])

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

      # make a diagnostic plot of cuts in the spatial direction
      self.spatialCuts()

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
      if visualize:
          # flat field has been medianed so most values should center on 1; want to be able to see finer-scale structure
          self.display.one(self.images['NormalizedFlat'], aspect='auto', vmin=0.95, vmax=1.05, scale='linear')
          self.display.run()

      plt.figure('normalized flat')
      ax = plt.subplot()
      ax.imshow(self.images['NormalizedFlat'].T, cmap='gray', aspect='auto', vmin=0.95, vmax=1.05)
      ax.set_xlabel('pixels')
      ax.set_ylabel('pixels')

      plt.savefig(normfilename)
      self.speak("saved normalized flat to {0}".format(normfilename))

      np.save(calibfilename, self.images)
      self.speak("saved normalized flat to calibration stamp {0}".format(calibfilename))

  def spatialCuts(self):

        ''' plot some cuts in the spatial direction for reference '''

        filename = os.path.join(self.directory, 'diagnosticCuts_{0}.pdf'.format(self.name))

        widths = self.trace.extractionwidths

        plt.figure('diagnistic cuts', figsize=(14,15), dpi=50)
        gs = plt.matplotlib.gridspec.GridSpec(5, 2, left=0.06, right=.98, bottom=0.04, top=0.98, hspace=0.17, wspace=0.15)
        cuts = {}
        cuts['stamp'] = plt.subplot(gs[:,0])

        image = self.images['science']/self.images['NormalizedFlat']
        extent = [self.saxis.min(), self.saxis.max(), self.waxis.min(), self.waxis.max()]

        cuts['stamp'].imshow(image, 
                             cmap='gray',
                             extent=extent,
                             interpolation='nearest',
                             aspect='auto',
                             zorder=0,
                             origin='lower',
                             vmin=0, vmax=np.percentile(image, 97))

        cuts['stamp'].set_xlabel('Spatial Direction (Pixels)')
        cuts['stamp'].set_ylabel('Dispersion Direction (Pixels)')

        # get the pixel indices of 5 rows in the cross-dispersion direction (don't use the first or last pixel rows though)
        take5 = np.linspace(0, len(self.waxis), 7)[1:-1]
        pixelcutinds = (np.round(take5)).astype(int)[::-1]

        axis, order = 0, 2
        for p, pix in enumerate(pixelcutinds):

            cuts[pix] = plt.subplot(gs[p, 1])

            profile = image.take(pix, axis=axis)
            cuts[pix].plot(self.saxis, profile, color='k', lw=2, alpha=.8)

            for width in widths:

                skymask = self.trace.skymask(width)
                extractmask = self.trace.extractionmask(width)

                sky = craftroom.twod.polyInterpolate(image, skymask == 0, order=order, visualize=False)

                cuts['stamp'].imshow(extractmask, 
                                     cmap=one2another('aquamarine', 'aquamarine', alphabottom=0.0, alphatop=1.0), 
                                     alpha=0.35, aspect='auto', extent=extent, origin='lower')
                cuts['stamp'].imshow(skymask, 
                                     cmap=one2another('deepskyblue', 'deepskyblue', alphabottom=0.0, alphatop=1.0),
                                     alpha=0.3, aspect='auto', extent=extent, origin='lower')

                notSky = skymask == 0
                extract = extractmask == 1
                
                notSkydata = notSky.take(pix, axis=axis)
                isExtractdata = extract.take(pix, axis=axis)
                extractdata = np.where(isExtractdata)[0]

                ok = notSkydata == 0
                med = np.median(profile[ok])
                mad = np.median(np.abs(profile[ok] - med))
                ok = ok*(np.abs(profile - med) <= 4*1.48*mad)

                fit = np.polynomial.polynomial.polyfit(self.saxis[ok], profile[ok], order)
                poly = np.polynomial.polynomial.Polynomial(fit)

                mincounts, maxcounts = np.min(profile), np.max(profile)

                # only need to plot the profile once per pixel cut; it doesn't change for different widths

                cuts[pix].plot(self.saxis, poly(self.saxis), lw=3, color='deepskyblue', alpha=0.3)
                cuts[pix].fill_between(self.saxis, 0, maxcounts, where=(notSkydata==0), color='deepskyblue', alpha=0.3)
                cuts[pix].fill_between(self.saxis, 0, maxcounts, where=(isExtractdata), color='aquamarine', alpha=0.3)

            cuts['stamp'].axhline(self.waxis[pix], self.saxis[0], self.saxis[-1], color='white', ls='--', lw=1.5, alpha=.5)

            cuts[pix].set_title('Cut on Pixel {0}'.format(int(self.waxis[pix])))
            cuts[pix].set_ylim(0,  np.percentile(profile, 96))
            cuts[pix].set_ylabel('Counts')
            cuts[pix].tick_params(labelbottom=False)

        cuts[pix].tick_params(labelbottom=True)
        cuts[pix].set_xlabel('Spatial Direction (Pixels)')

        plt.savefig(filename)
        self.speak("saved plot of diagnostic cuts in spatial direction")
        plt.close()
        

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
    self.extracted = np.load(self.extractedFilename, allow_pickle=True)[()]
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

        #self.extracted['original'] = self.intermediates['original']



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

            # want to save some of the intermediates stuff into extracted for later, just in case
            self.extracted[width]['intermediates'] = {}
            self.extracted[width]['intermediates']['skyMask'] = self.intermediates[width]['skyMask']
            self.extracted[width]['intermediates']['extractMask'] = self.intermediates[width]['extractMask']

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

                # pick out the exposure prefixes that will get plotted as a diagnostic
                fileprefixes = self.obs.fileprefixes['science']
                indices = np.floor(np.linspace(0, len(fileprefixes)-1, 6)).astype('int')

                self.speak('using FWHM modification')
                mask = (self.images['badpixels']-1)*-1
                subimage = (image/self.images['NormalizedFlat']*mask)#[:,boxcuts[0]:boxcuts[-1]]
                disp, crossdisp = subimage.shape
                FWHM = np.zeros(disp)
                edge1, edge2 = np.zeros(disp), np.zeros(disp)

                subarray = self.instrument.extractiondefaults['spatialsubarray']
                lastr1, lastr2 = subarray, subarray+width
                for i in range(disp):

                    # make a psf by cross-dispersion column; multiply by extraction mask so neighboring spectra are not included
                    psf = subimage[i]*self.intermediates[width]['extractMask'][i]
                    # create a spline to fit to the psf; this will give the FWHM roots
                    spline = scipy.interpolate.UnivariateSpline(range(len(psf)), psf-np.max(psf)/2.)
                    #print(i)
                    #print(psf)
                    #ipyprint(spline.roots())
                    try:
                        r1, r2 = spline.roots() # find the roots
                        lastr1, lastr2 = r1, r2
                    except(ValueError):
                        #print(i, width, 'roots: ', spline.roots(), 'using: ', lastr1, lastr2)
                        r1, r2 = lastr1, lastr2

                    FWHM[i] = r2-r1
                    edge1[i] = r1
                    edge2[i] = r2

                # up to here edges will be basically the same for every extraction width;
                # extend edges out by width/2 times FWHM/2
                edge1 -= (FWHM/2.) * (width/2.)
                edge2 += (FWHM/2.) * (width/2.)

                # smooth edges
                edge1spline = scipy.interpolate.UnivariateSpline(range(len(edge1)), edge1)
                edge1smooth = edge1spline(range(len(edge1)))
                edge2spline = scipy.interpolate.UnivariateSpline(range(len(edge2)), edge2)
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
                self.intermediates[width]['edge1smooth'] = edge1smooth
                self.intermediates[width]['edge2smooth'] = edge2smooth

                self.extracted[width]['median_width'] = np.median(edge2smooth - edge1smooth)

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

                # comment this next line out if you don't need it; it makes the extraction a lot slower
                #self.extracted[width]['intermediates']['subtracted'] = self.intermediates[width]['subtracted']

                # this is a plotting tool to go along with the FWHM modification
                #if self.exposureprefix in [fileprefixes[i] for i in indices]:
                #    diagnosticFilename = os.path.join(self.directory, 'extracted_diagnostic_{0}_{1}px.pdf'.format(self.exposureprefix, width))
                #    plt.figure('diagnostic')
                #    plt.imshow(image/self.images['NormalizedFlat' ] - self.intermediates[width]['sky'], aspect='auto', interpolation='none', vmin=0, vmax=50000)
                #    plt.plot(edge1smooth, range(len(self.waxis)), color='C1')
                #    plt.plot(edge2smooth, range(len(self.waxis)), color='C1')
                #    plt.colorbar()
                #    plt.savefig(diagnosticFilename)
                #    plt.close('diagnostic')

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

                # add intermediates stuff to extracted for later inspection
                #self.extracted[width]['intermediates']['sky'] = self.intermediates[width]['sky']


                # KLUDGE -- am I definitely catching this error later?
                #for k in ['centroid', 'width', 'peak']:
                #    assert(np.isfinite(self.extracted[width][k]).all())


            else:
                self.extracted[width]['raw_counts'] = np.nansum(self.intermediates[width]['extractMask']*image/self.images['NormalizedFlat'], self.sindex)
                # this is a kludge, to make the plotting look better for arcs
                self.intermediates[width]['sky']  = np.zeros_like(self.intermediates['original'])# + np.percentile(self.extracted[width]['raw_counts'] , 1)

            # diagnostic: saves a 
            #import pickle
            #takefive = int(len(self.obs.fileprefixes['science'])/5)
            #if self.exposureprefix in self.obs.fileprefixes['science'][::takefive]:# & (width == 6.0):
            #    pickle.dump(self.intermediates, open('/home/hdiamond/LHS1140/from_extraction/intermediates_'+self.obs.night.name+'_'+self.name+'_'+self.exposureprefix+'_'+str(width)+'px.p', 'wb'))

            #import sys
            #sys.exit("Breaking here. Check it out.")

        if arc==False:
            self.visualizeExtraction()


        # save the extracted spectra for all apertures
        np.save(self.extractedFilename, self.extracted)
        self.speak('saved extracted spectra to {0}'.format(self.extractedFilename))

    # create a supersampled version of the extracted spectra
    #if False:#arc == False:
    #    self.interpolate(remake=remake)

    # return the extracted spectrum
    return self.extracted
 
  def setupVisualization(self):
        self.thingstoplot = ['raw_counts']#['raw_counts_optext', 'sky', 'width',  'raw_counts']
        height_ratios = np.ones(len(self.thingstoplot) + 2)
        suptitletext = '{}, {}'.format(self.name, self.exposureprefix)
        try:
            self.suptitle.set_text(suptitletext)
            self.plotted
        except AttributeError:
            self.figure = plt.figure(figsize=(20,12), dpi=50)
            gs = plt.matplotlib.gridspec.GridSpec(len(self.thingstoplot) + 2, self.trace.numberofapertures,
                                                    hspace=0.08, wspace=0.05,
                                                    height_ratios = height_ratios)
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
            
                # row of 2D plots showing extraction stamp with extraction mask and sky mask for each width
                self.ax[width]['apertures'] = plt.subplot(gs[-2, i], sharex=sharex, sharey=shareapertures)
                sharex = self.ax[width]['apertures']
                shareapertures = self.ax[width]['apertures']
                self.ax[width]['apertures'].tick_params(labelbottom=False)

                # row of 2D plots showing extraction stamp with extraction mask and sky mask with the sky background subtracted for each width
                self.ax[width]['subtracted'] = plt.subplot(gs[-1, i], sharex=sharex, sharey=sharesubtracted)
                sharesubtracted = self.ax[width]['subtracted']
                self.ax[width]['subtracted'].set_xlabel('Dispersion Direction (Pixels)')
                self.ax[width]['apertures'].set_ylim(*craftroom.oned.minmax(self.saxis))

                # set the y labels
                if i == 0:
                    self.ax[width]['apertures'].set_ylabel('extraction apertures')
                    self.ax[width]['subtracted'].set_ylabel('sky-subtracted, and coarsely rectified')
                else:
                    self.ax[width]['subtracted'].tick_params(labelleft=False)
                    self.ax[width]['apertures'].tick_params(labelleft=False)


                for j, thing in enumerate(self.thingstoplot):
                    self.ax[width][thing] = plt.subplot(gs[j,i], sharex=sharex, sharey=sharey[thing])
                    self.ax[width][thing].locator_params(nbins=3)
                    sharey[thing] = self.ax[width][thing]

                    if j == 0:
                        self.ax[width][thing].set_title('{} pix.'.format(width))
                    if i == 0:
                        self.ax[width][thing].set_ylabel(thing)
                    else:
                        self.ax[width][thing].tick_params(labelleft=False)
                    self.ax[width][thing].tick_params(labelbottom=False)

                    if thing == 'width':
                        self.ax[width][thing].set_ylim(0,np.max(self.trace.extractionwidths)/3.5)
                    if thing == 'centroid':
                        self.ax[width][thing].set_ylim(np.min(self.trace.traceCenter(self.waxis))-5, np.max(self.trace.traceCenter(self.waxis))+5)
                    if thing == 'raw_counts':
                        self.ax[width][thing].set_ylim(0, np.percentile(self.extracted[width]['raw_counts'], 99)*1.5)
                    if thing == 'raw_counts_optext':
                        self.ax[width][thing].set_ylim(0, np.percentile(self.extracted[width]['raw_counts_optext'], 99)*1.5)

  def visualizeExtraction(self):
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

  def plotExtracted(self, width, ax=None):
        '''plot the extract spectra and parameters'''

        if ax is None:
            ax = self.ax[width]

        for j, thing in enumerate(self.thingstoplot):
            if thing in self.extracted[width].keys():
                try:
                    self.plotted[width][thing][0].set_data(self.extracted['w'], self.extracted[width][thing])
                except KeyError:
                    self.plotted[width][thing] = \
                        ax[thing].plot(
                            self.extracted['w'],
                            self.extracted[width][thing],
                            linewidth=2, color='k', alpha=0.8)

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
          self.extracted = np.load(self.extractedFilename, allow_pickle=True)[()]

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
        '''
        Interpolate the spectra onto a common (uniform) wavelength scale.
        This includes both the raw_counts, as well as all the auxiliary
        wavelength-axis arrays (sky, peak, centroid, width) too. It saves
        these resampled spectra to their own file.

        Parameters
        ----------
        remake : bool
            Should we remake the resampled arrays,
            even if some saved files already exist?

        shift : bool
            Should we apply a shift/stretch to the wavelengths
            when doing the resampling? A general path will be
            to extract first with shift=False, then cross-correlate
            multiple spectra (possibly in both times and stars)
            to estimate time-dependent stretches, then rerun
            this function with shift=True, so the resampling
            uses the updated wavelengths.

        '''


        # these are things where we care about the sum matching between extracted and supersampled
        self.additivekeys = ['raw_counts', 'sky']

        # these are things where we want the individual values matching between extracted and supersampled
        self.intrinsickeys = ['centroid', 'width', 'peak']

        # these are all the keys that will be supersampled
        self.keys = self.additivekeys + self.intrinsickeys


        try:
            # try just loading the supersampled spectra from files
            assert(remake == False)
            self.supersampled = np.load(self.supersampledFilename, allow_pickle=True)
            self.speak('loaded supersampled spectrum from {0}'.format(self.supersampledFilename))
        except (AssertionError, IOError):
            # otherwise, make some new supersampled spectra
            if shift:

                try:
                    # do we already know how to stretch spectra in this aperture?
                    self.stretches
                except AttributeError:
                    # load the spectral stretches that were estimated via a WavelengthRecalibrator
                    stretchfilename = os.path.join(self.mask.reducer.extractionDirectory,  'spectralstretch.npy')
                    self.stretches = np.load(stretchfilename, allow_pickle=True)[()]
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

                # the wavelength originally estimated for the pixels
                originalwavelength = self.extracted['wavelength']

                # pull out the ingredients for stretching/shifting the wavelengths
                star = self.name
                midpoint = self.stretches['midpoint']
                coefficients = self.stretches['stretch'][star][self.exposureprefix], self.stretches['shift'][star][self.exposureprefix]
                # calculate dw = originalwavelength - newwavelength  = f(originalwavelength - midpoint)
                dw = np.polyval(coefficients, originalwavelength - midpoint)
                wavelength = originalwavelength - dw

                # quote by how much the wavelengths were nudged
                phrase = 'dw = {:.4}x(w - {midpoint}){:+.4}'.format(*coefficients, midpoint=midpoint)
                self.speak('nudge wavelengths for {} by {}'.format(self.exposureprefix, phrase))
            else:
                wavelength = self.extracted['wavelength']

            # what's the original pixel number for each extracted bin?
            pixelnumber = self.extracted['w']

            # set up a (generally higher resolution) common wavelength grid onto which everything will be interpolated
            try:
                # does the supersampled grid already exist?
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

                # store our shared wavelength grid
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
                    ysupersampled[yclosetonan] = np.nan

                    # decide whether we are summing or averaging this array
                    if key in self.additivekeys:
                        # fluxconservingresample by default gives a sum
                        self.supersampled[combinedkey] = ysupersampled
                    elif key in self.intrinsickeys:
                        # divide by the number of pixels to get down to an average
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
