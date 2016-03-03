from imports import *
from Tools import *
colors = dict(He='lightsalmon', Ne='red', Ar='deepskyblue')

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
    self.display = self.calib.display#zachods9('aperture', xsize=800, ysize=200, rotate=90)
    self.setup(x,y)
    self.createCalibStamps()
    self.createTrace()
    self.createSkyApertures()
    self.createWavelengthCal()

  def stampFilename(n):
	'''Spit out the right stamp filename for this aperture, cut out of CCDn.'''
	return self.directory + 'stamp{0:04}.fits'.format(n)

  def setup(self,x,y):
    '''Setup the basic geometry of the aperture.'''
    if self.obs.instrument == 'LDSS3C':
      self.x = x
      self.y = y
      self.maskWidth = (self.obs.skyWidth + self.obs.skyGap)*2 #+20 +
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
      self.w = self.y_sub - self.ystart
      self.s = self.x_sub - self.xstart
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
    filename = self.directory + 'trace_{0}.npy'.format(self.name)
    try:
      traceCoeff, width = np.load(filename)
      self.speak("loaded previously fitted trace from {0}".format(filename))
    except IOError:
      self.speak("fitting master science image to determine trace parameters")
      converged = False

      # start with a width guess and assuming the star is in middle of aperture
      old_width = self.obs.widthGuess
      traceCoeff = [np.max(self.saxis)/2]

      # iterate through, trying to improve
      width = old_width
      attempts = 0
      while(converged == False):

          self.speak('on iteration {0} of trace-fitting'.format(attempts))

          # set the trace parameters right now
          self.traceCenter = np.polynomial.polynomial.Polynomial(traceCoeff)
          self.traceWidth = width

          self.speak(" trying trace coefficients of {0}, width of {1:.1f}".format(traceCoeff, width))

          # consider points within 2.5 sigma of trace as the star

          distancefromtrace = np.abs((self.s - self.traceCenter(self.w)))
          buffer = 2.5*self.traceWidth
          considerstar = distancefromtrace < buffer
          sortofsky = (distancefromtrace < buffer + self.obs.skyWidth)
          considersky = (considerstar == False)*sortofsky

          # estimate a rough 1D spectrum
          flattened = self.images['Science']/self.images['NormalizedFlat']
          roughSky1d=np.average(flattened,
                        axis=self.sindex,
                        weights=considersky)
          reshaped = roughSky1d.reshape((self.waxis.shape[0],1))
          self.images['Sky'] = reshaped*np.ones_like(flattened)
          self.images['Subtracted'] = flattened - self.images['Sky']

          # use the roughly subtracted image to fit centroids
          fluxWeightedCentroids = np.average(self.s,
                                axis=self.sindex,
                                weights=considerstar*self.images['Subtracted'])

          fluxWeightedWidths = np.sqrt(np.average(self.s**2,
                                axis=self.sindex,
                                weights=considerstar*self.images['Subtracted']) - fluxWeightedCentroids**2)

          fluxWeightedWidths[np.isfinite(fluxWeightedWidths) == False] = np.inf

          # fit a polynomial to the ridge of the spectrum
          traceCoeff = np.polynomial.polynomial.polyfit(self.waxis, fluxWeightedCentroids, self.obs.traceOrder, w=1.0/fluxWeightedWidths**2)

          fit_width = np.median(fluxWeightedWidths)
          width = np.minimum(fit_width, self.obs.widthGuess)

          # create a rough LSF
          self.images['RoughLSF'] = np.exp(-0.5*((self.s - self.traceCenter(self.w))/self.traceWidth)**2)

          converged = np.abs(fit_width - old_width) < 0.01
          self.speak( "  {0} -> {1}, converged = {2}".format(old_width, width, converged))
          old_width = width



          self.display.rgb(self.images['Subtracted'], self.images['RoughLSF'], considersky)
          self.display.one(self.images['Science'], clobber=False)
          self.display.tile('rows')
          #self.input('continue?')
          attempts += 1

          if attempts > 10:
            old_width = int(self.input('what width should we try?'))
            width = old_width
            attempts=0

      assert("n" not in self.input("like trace?").lower())
      np.save(filename, (traceCoeff, width))
      self.speak("saved trace parameters to {0}".format( filename))

    self.traceCenter = np.polynomial.polynomial.Polynomial(traceCoeff)
    self.traceWidth = width
    self.images['RoughLSF'] = np.exp(-0.5*((self.s - self.traceCenter(self.w))/self.traceWidth)**2)


  def createSkyApertures(self, visualize=True):
    '''Let user select the best sky apertures.'''
    self.speak("setting up the sky apertures")
    filename = self.directory + 'skyMask_{0}.npy'.format(self.name)
    try:
      self.images['skyMask'] = np.load(filename)
      self.speak("loaded sky apertures from {0}".format(filename))
    except:
      finished = False
      plt.figure('sky apertures', figsize=(10,10), dpi=100)
      i = zachopy.iplot.iplot(3,10)
      self.aximage = i.subplot(1,0,rowspan=2,colspan=8)
      self.axskyspectrum = i.subplot(0,0,colspan=8,sharex=self.aximage)
      self.axskyprofile = i.subplot(1,8,rowspan=2,colspan=2)
      axes = [self.aximage, self.axskyspectrum, self.axskyprofile]
      mask = np.zeros_like(self.images['Science'])
      extent=[self.waxis.min(), self.waxis.max(), self.saxis.min(), self.saxis.max()]
      finishedplotting = True
      first = True
      while(finished == False):


        if first == False:
          # have user select a sky region
          self.speak(":-D please click the start of a sky region")
          clicks = i.getMouseClicks(n=2)

          # clear the axes
          for a in axes:
            a.cla()


        # display the image
        self.aximage.imshow(np.transpose(np.log(self.images['Science'])), cmap='gray', \
                      extent=extent, \
                      interpolation='nearest', aspect='auto', \
                      vmin=np.log(1), vmax=np.log(np.nanmax(self.images['Science'])*0.01))
        self.aximage.imshow(np.transpose(mask), alpha=0.1, cmap='winter_r', \
                    extent=extent, \
                    interpolation='nearest', aspect='auto')

        # overlay the trace
        self.aximage.plot(self.waxis, self.traceCenter(self.waxis), color='blue', alpha=0.3, linewidth=4)
        self.aximage.set_xlim(self.waxis.min(), self.waxis.max())


        if first == False:
          # calculate offsets from the trace
          offsets = (clicks[0].ydata - self.traceCenter(clicks[1].xdata), clicks[1].ydata - self.traceCenter(clicks[1].xdata))
          bottom = np.min(offsets)
          top = np.max(offsets)

          # display the most recent
          self.aximage.plot(self.waxis, self.traceCenter(self.waxis) + bottom , color='green', alpha=0.3, linewidth=4)
          self.aximage.plot(self.waxis, self.traceCenter(self.waxis) + top , color='green', alpha=0.3, linewidth=4)
          mask[(self.s > self.traceCenter(self.w) + bottom) * (self.s < self.traceCenter(self.w) + top)] += 1.0
          mask = mask > 0

          ma = np.ma.MaskedArray(self.images['Science'], mask==0)
          skyspectrum = np.ma.median(ma, self.sindex)
          self.axskyspectrum.plot(skyspectrum)
          click = clicks[-1]
          self.axskyprofile.cla()
          self.axskyprofile.plot(self.images['Science'][click.xdata,:], self.saxis)
          self.axskyprofile.plot((self.images['Science']*mask)[click.xdata,:], self.saxis, linewidth=3)
          self.axskyprofile.set_xlim((self.images['Science']*mask)[click.xdata,:].min(), (self.images['Science']*mask)[click.xdata,:].max()*2)
          self.axskyprofile.set_ylim(self.saxis.min(), self.saxis.max())
          plt.draw()
          self.speak("Are you happy with the sky subtraction apertures? (default = no)")
          answer = self.input("  (y)es, (n)o, (r)edo")
          if "y" in answer:
            finished = True
          elif "r" in answer:
            mask *= 0
          else:
            finished = False
        first = False
      self.images['skyMask'] = mask
      np.save(filename, self.images['skyMask'])
      self.speak("saved a sky mask to {0}".format(filename))

      # cross correlate my arc spectra with reference, to find rough offset
  def findRoughShift(self, wavelength_ids, visualize=True):
    '''using a list of pixels matched to wavelengths, find the rough
            offset of the this arc relative to the standard slits/setting'''

    self.speak("cross correlating arcs with known wavelengths")

    # create a plot showing how well the lines match
    if self.visualize:
        figure_waverough = plt.figure('wavelength rough offset',
                                            figsize=(10,5), dpi=100)
        gs = plt.matplotlib.gridspec.GridSpec(  3,2,
                                                bottom=0.15, top=0.85,
                                                hspace=0.1,wspace=0,
                                                width_ratios=[1, .5])
        ax_waverough, ax_wavecor = {}, {}
        sharer, sharec = None, None
        for i, e in enumerate(self.elements):
            ax_waverough[e] = plt.subplot(gs[i,0],sharex=sharer)
            ax_wavecor[e] = plt.subplot(gs[i,1], sharex=sharec)
            sharer, sharec = ax_waverough[e], ax_wavecor[e]

    correlationfunctions = {}
    for count, element in enumerate(self.elements):
      # pull out the flux and pixel coordinates from extract arc spectrum
      flux = self.arcs[element]['raw_counts']
      x = self.arcs[element]['w']

      # find peaks in this spectrum
      xPeak, yPeak = zachopy.oned.peaks(self.waxis, flux)

      # pull out the line identifications that match this element
      this = []
      for i in range(len(wavelength_ids)):
        if element in wavelength_ids['name'][i]:
          this.append(i)
      self.obs.binning = 2
      # create fake spectra using the line positions (both reference + new)
      myPeaks, theirPeaks = np.zeros(len(x)), np.zeros(len(x))
      for i in range(len(this)):
        center = wavelength_ids['pixels'][this[i]]/self.obs.binning
        blob = 2
        theirPeaks += np.exp(-0.5*((x-center)/blob)**2)
      for i in range(len(xPeak)):
        center = xPeak[i]
        blob = 2
        myPeaks += np.exp(-0.5*((x-center)/blob)**2)*np.log(yPeak[i])

      correlationfunctions[element] = np.correlate(myPeaks, theirPeaks, 'full')
      if self.visualize:

        ax_waverough[element].plot(x, myPeaks/myPeaks.max(),
                                        label='extracted', alpha=0.5,
                                        color=colors[element])
        ax_waverough[element].set_ylim(0, 2)
        normcor = correlationfunctions[element]
        normcor /= np.max(correlationfunctions[element])
        ax_wavecor[element].plot(normcor,
                                    label=element, alpha=0.5,
                                    color=colors[element])
        for a in [ax_wavecor, ax_waverough]:
            plt.setp(a[element].get_xticklabels(), visible=False)
            plt.setp(a[element].get_yticklabels(), visible=False)

      if count == 0:
        corre =  correlationfunctions[element]
      else:
        corre = corre*correlationfunctions[element]


    # find the peak of the combined correlation function
    peakoffset = np.where(corre == corre.max())[0][0] - len(x)
    # to convert: len(x) - xPeak = x + peakoffset


    if self.visualize:
      for element in self.elements:
        for i in range(len(wavelength_ids)):
          if element in wavelength_ids['name'][i]:
            center = peakoffset  + wavelength_ids['pixels'][i]/self.obs.binning
            ax_waverough[element].axvline(center ,
                                          alpha=0.25, color='black')
        ax_wavecor[element].plot(corre/np.max(corre),
                            label='combined', alpha=0.25, color='black')
        ax_waverough[element].set_ylabel(element)
      ax = ax_waverough[self.elements[-1]]
      plt.setp(ax.get_xticklabels(), visible=True)
      ax.set_xlabel('Pixel Position')
      ax_wavecor[self.elements[0]].set_title('cross correlation peaks at \n{0} pixels ({1}x{1} binned pixels)'.format(peakoffset, self.obs.binning))

      ax_waverough[self.elements[0]].set_title('Coarse Wavelength Alignment\nfor (%0.1f,%0.1f)' % (self.x, self.y))

      figure_waverough.savefig(self.directory + 'roughWavelengthAlignment_{0}.pdf'.format(self.name))


    return peakoffset

  def createWavelengthCal(self, wavelengthIDFile=None, visualize=True):
    '''Populate the wavelength calibration for this aperture.'''
    self.speak("populating wavelength calibration")

    self.elements = ['He', 'Ne','Ar']
    # if the wavelength cal already exists, simply reload it
    filename = self.directory + 'waveCal_{0}.npy'.format(self.name)
    try:
      # load the the thing!
      self.waveCalCoef = np.load(filename)
    except:
      self.speak("estimating wavelength calibration from extracted arc spectra")
      # extract the arc lamp spectra

      # extract a spectrum from the master image for each lamp
      self.arcs = {}
      for element in self.elements:
        self.arcs[element] = self.extract(image=self.images[element], arc=True)

      # load wavelength identifications (matched to pixels)
      if wavelengthIDFile is None:
          wavelengthIDFile = self.obs.wavelength2pixelsFile
      wavelength_ids = astropy.io.ascii.read(wavelengthIDFile)

      # load the complete list of wavelengths
      allwavelengths = astropy.io.ascii.read(self.obs.wavelengthsFile)

      # perform a first rough alignment, to make pixels
      peakoffset = self.findRoughShift(wavelength_ids)

      notconverged = True
      #while(notconverged):
      if True:
          # empty lists, which we'll fill with pixels matched to wavelengths
          pixel, wavelength, emissioncolors = [],[], []

          wavelengthorder = 3

          # create a temporary calibration to match reference wavelengths to reference pixels (so we can include additional wavelengths not recorded in the dispersion solution file)

          p = np.polynomial.polynomial
          coef = p.polyfit(
                (wavelength_ids['wave2']),
                wavelength_ids['pixels']/self.obs.binning +  peakoffset,
                wavelengthorder)# + wavelength_ids['wave2'])/2.0
          theirWavelengthstoMyPixels = p.Polynomial(coef)


          # treat the arc lamps separately
          for count, element in enumerate(self.elements):

            # the pixel spectrum self.extracted from this arc lamp
            flux = self.arcs[element]['raw_counts']
            x = self.arcs[element]['w']

            # identify my peaks
            xPeak, yPeak = zachopy.oned.peaks(self.waxis, flux)


            # pull out the wavelengths from the complete file
            theirWavelengths = []
            for i in range(len(allwavelengths)):

              if element in allwavelengths['name'][i]:
                theirWavelengths.append(allwavelengths['wavelength'][i])
            theirWavelengths = np.array(theirWavelengths)
            print element
            print theirWavelengths
            # put those peaks onto my pixel scale
            theirPeaksOnMyPixels = theirWavelengthstoMyPixels(theirWavelengths)

            iMine, iTheirs =[],[]

            # pull out my peaks
            myPeaksOnMyPixels = xPeak
            for i in range(len(theirPeaksOnMyPixels)):
              # find my closest peak to theirs
              distance = xPeak - theirPeaksOnMyPixels[i]
              closestdistance = (np.abs(distance)).min()

              # if the peak matches closely, use it
              if closestdistance < 10:
                iMine.append(np.where(np.abs(distance) == closestdistance)[0][0])
                iTheirs.append(i)

            initialCoeff = np.polynomial.polynomial.polyfit(xPeak[iMine], theirWavelengths[iTheirs], 2)
            initialFit = np.polynomial.polynomial.Polynomial(initialCoeff)
            pixel.extend(xPeak[iMine])
            wavelength.extend(theirWavelengths[iTheirs])
            emissioncolors.extend([colors[element]]*len(iMine))
            count += 1

          fit = np.polynomial.polynomial.Polynomial(np.polynomial.polynomial.polyfit(pixel, wavelength, wavelengthorder) )
          residual = wavelength - fit(pixel)
          good = np.abs(residual) < np.median(np.abs(residual))*3
          bad = ~good
          self.waveCalCoef  = np.polynomial.polynomial.polyfit(pixel, wavelength, wavelengthorder, w=good)
          fit = np.polynomial.polynomial.Polynomial(self.waveCalCoef )
          pixel, wavelength = np.array(pixel), np.array(wavelength)

          if self.visualize:
            # plot to make sure the wavelength calibration makes sense
            figure_wavelengthcal = plt.figure('wavelength calibration')
            #gs = plt.matplotlib.gridspec.GridSpec(

            interactivewave = zachopy.iplot.iplot(4,1,
                    height_ratios=[.2, .2, .4, .2], hspace=0)


            # does the wavelength2pixel code work
            ax_w2p = interactivewave.subplot(0)

            # do the lamp spectra overlap?
            ax_walign = interactivewave.subplot(1, sharex=ax_w2p)

            # what is the actual wavelength calibration
            ax_wcal = interactivewave.subplot(2, sharex=ax_w2p)

            # what are the residuals from the fit
            ax_wres = interactivewave.subplot(3, sharex=ax_w2p)

            for ax in [ax_w2p, ax_walign, ax_wcal]:
                plt.setp(ax.get_xticklabels(), visible=False)
            ax_w2p.set_title('Wavelength Calib. for Aperture (%0.1f,%0.1f)' % (self.x, self.y))

            kw = dict(color=emissioncolors, marker='o')

            # plot the backwards calibration
            ax_w2p.scatter(wavelength_ids['pixels']/self.obs.binning + peakoffset, wavelength_ids['wave2'], **kw)
            xfine = np.linspace(min(wavelength_ids['wave2']), max(wavelength_ids['wave2']))
            ax_w2p.plot(theirWavelengthstoMyPixels(xfine), xfine)

            # plot the overlap of the lamp spectra
            for element in self.elements:

                ax_walign.plot(self.waxis, self.arcs[element]['raw_counts'],
                                    color=colors[element], alpha=0.5)
                ax_walign.set_yscale('log')


            for i,tw in enumerate(allwavelengths['wavelength']):
                name = allwavelengths['name'][i][0:2]
                if name in self.elements:
                    ax_walign.axvline(theirWavelengthstoMyPixels(tw), ymin=0.9,
                                    color=colors[name],
                                    alpha=0.5)

            # plot the new calibration
            ax_wcal.scatter(pixel, wavelength, **kw)
            ax_wcal.plot(pixel, fit(pixel), alpha=0.5, color='black')
            ax_wcal.set_ylabel('Wavelength (angstroms)')

            ax_wres.set_ylabel('Residuals')
            ax_wres.scatter(pixel[good], wavelength[good] - fit(pixel[good]), **kw)
            kw['marker'] = 'x'
            ax_wres.scatter(pixel[bad], wavelength[bad] - fit(pixel[bad]), **kw)
            ax_wres.set_xlabel('Pixel # (by python rules)')
            ax_wres.set_xlim(min(pixel), max(pixel))


            if 'n' in self.input('are you okay with wavelength cal?').lower():
                self.createWavelengthCal('/Users/zkbt/Dropbox/code/mosasaurus/data/vph-red_wavelength_identifications_blueslit.txt')

            #clicks = interactivewave.getMouseClicks(n=2)
            #print clicks

            figure_wavelengthcal.savefig(self.directory + 'wavelengthCalibration_{0}.pdf'.format(self.name))

      np.save(filename, self.waveCalCoef)
      self.speak("saved wavelength calibration to {0}".format( filename))
    self.wavelengthCalibrate = np.polynomial.polynomial.Polynomial(self.waveCalCoef )

  def ones(self):
    '''Create a blank array of ones to fill this aperture.'''
    return np.ones_like(self.images['Science'])

  def extract(self, n=0, image=None, subtractsky=True, arc=False, cosmics=False, remake=False):
    '''Extract the spectrum from this aperture.'''

    self.extractedFilename = self.directory + 'extracted{0:04}.npy'.format(n)
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
          self.extracted['wavelength'] = self.wavelengthCalibrate(self.waxis)
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
