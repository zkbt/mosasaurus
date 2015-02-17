from imports import *

class Aperture(Talker):
  '''Aperture objects store all information about individual apertures (e.g. slits).'''
  def __init__(self, x, y, mask, **kwargs):
    '''Initialize an aperture, taking a position on the detector as input.'''

    # decide whether or not this Reducer is chatty
    Talker.__init__(self, **kwargs)
    self.visualize = True
    self.mask = mask
    self.calib = self.mask.calib
    self.obs = self.calib.obs
    self.display =  zachopy.display.ds9('aperture')
    self.setup(x,y)
    self.createCalibStamps()
    self.createTrace()
    self.createSkyApertures()
    self.createWavelengthCal()


  def setup(self,x,y):
    '''Setup the basic geometry of the aperture.'''
    if self.obs.instrument == 'LDSS3C':
      self.x = x
      self.y = y
      self.maskWidth = (self.obs.skyWidth +20 + self.obs.skyGap)*2
      self.ystart = np.maximum(y - self.obs.blueward, 0)
      self.yend = np.minimum(y + self.obs.redward, self.obs.ysize)
      self.xstart = np.maximum(x - self.maskWidth, 0)
      self.xend = np.minimum(x + self.maskWidth, self.obs.xsize)
      # remember python indexes arrays by [row,column], which is opposite [x,y]
      x_fullframe, y_fullframe = np.meshgrid(np.arange(self.calib.images['Science'].shape[1]),np.arange(self.calib.images['Science'].shape[0]))
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
      interesting = ['Science' ,'WideFlat', 'He', 'Ne', 'Ar', 'BadPixels', 'Dark', 'Bias']
      for k in interesting:
        self.images[k] = self.stamp(self.calib.images[k])

      self.speak('these are the stamps before interpolating over bad pixels')
      self.displayStamps(self.images, keys = ['Science', 'WideFlat', 'BadPixels'])
      self.input('', prompt='(press return to continue)')
      for k in interesting:
          if k != 'BadPixels':
              self.images[k] = zachopy.twod.interpolateOverBadPixels(self.images[k], self.images['BadPixels'])

      self.speak('and these are they after interpolating over bad pixels')
      self.displayStamps(self.images, keys = ['Science', 'WideFlat', 'BadPixels'])
      self.input('', prompt='(press return to continue)')

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


      # create a sky subtracted science image, for estimate shape of trace
      roughSky1d = np.median(self.images['Science']/self.images['NormalizedFlat'], self.sindex)
      #zachopy.twod.estimateBackground(self.images['Science'], axis=self.sindex)
      self.images['Sky'] = np.ones_like(self.images['Science'])*roughSky1d.reshape((self.waxis.shape[0],1))
      self.images['Subtracted'] = (self.images['Science'] - self.images['Sky']*self.images['NormalizedFlat'])#np.maximum(, 0)
        # (the edges of the slit might not be exposed -- they'd drop to way negative without the np.maximum statement


      if visualize:
        try:
            self.ax.cla()
        except:
            self.fi, self.ax = plt.subplots(1,1)
        self.ax.cla()
        self.ax.plot(self.waxis, envelope)
        self.ax.plot(self.waxis, spline(self.waxis))

      if visualize:
        self.displayStamps(self.images, keys = ['Science', 'Sky', 'Subtracted', 'WideFlat', 'NormalizedFlat', 'BadPixels'])
        assert("n" not in self.input('Do you like the calibration stamps?').lower())

      np.save(filename, self.images)
      self.speak("saved calibration stamps to {0}".format( filename))

  def displayStamps(self, images, keys=None):
    '''Display stamps relevant to this aperture in ds9.'''
    if keys is None:
        keys = images.keys()

    self.display.tile('row')
    for i in range(len(keys)):
      self.display.replace(np.transpose(images[keys[i]]), i)

  def createTrace(self):
    '''Fit for the position and width of the trace.'''
    self.speak("populating the trace parameters")
    filename = self.directory + 'trace_{0}.npy'.format(self.name)
    try:
      traceCoeff, width = np.load(filename)
      self.speak("loaded previously fitted trace from {0}".format(filename))
    except:
      self.speak("fitting to the master science image to determine the trace parameters")
      converged = False
      old_width = 10
      traceCoeff = [np.max(self.saxis)/2]
      width = old_width
      while(converged == False):

        self.traceCenter = np.polynomial.polynomial.Polynomial(traceCoeff)
        self.traceWidth = width
        weighting = np.abs((self.s - self.traceCenter(self.w))/self.traceWidth) < 2.5
        # estimates the centroids of the spectrum in the spatial direction
        fluxWeightedCentroids = (self.s*self.images['Subtracted']*weighting).sum(self.sindex)/(self.images['Subtracted']*weighting).sum(self.sindex)
        fluxWeightedWidths = np.sqrt((self.s**2*self.images['Subtracted']*weighting).sum(self.sindex)/(self.images['Subtracted']*weighting).sum(self.sindex) - fluxWeightedCentroids**2)
        fluxWeightedWidths[np.isfinite(fluxWeightedWidths) == False] = np.inf

        # fit a polynomial to the ridge of the spectrum
        traceCoeff = np.polynomial.polynomial.polyfit(self.waxis, fluxWeightedCentroids, self.obs.traceOrder, w=1.0/fluxWeightedWidths**2)
        fit_width = np.median(fluxWeightedWidths)
        width = np.minimum( fit_width, 10)



        self.images['RoughLSF'] = np.exp(-0.5*((self.s - self.traceCenter(self.w))/self.traceWidth)**2)
        converged = np.abs(fit_width - old_width) < 0.01
        self.speak( "{0} -> {1}, converged = {2}".format(old_width, width, converged))
        old_width = width



        self.display.rgb(self.images['Subtracted'], self.images['RoughLSF'], weighting)
        self.display.one(self.images['Science'], clobber=False)
      assert("n" not in self.input("like trace?").lower())
      np.save(filename, (traceCoeff, width))
      self.speak("saved trace parameters to {0}".format( filename))

    self.traceCenter = np.polynomial.polynomial.Polynomial(traceCoeff)
    self.traceWidth = width
    self.images['RoughLSF'] = np.exp(-0.5*((self.s - self.traceCenter(self.w))/self.traceWidth)**2)
    self.display.rgb(self.images['Subtracted'], self.images['RoughLSF'], self.images['Sky'])
    self.display.one(self.images['Science'], clobber=False)

  def createSkyApertures(self, visualize=True):
    '''Let user select the best sky apertures.'''
    print "     setting up the sky apertures"
    filename = self.directory + 'skyMask_{0}.npy'.format(self.name)
    try:
      self.images['skyMask'] = np.load(filename)
      print "         loaded sky apertures from {0}".format(filename)
    except:
      finished = False
      plt.figure('sky apertures')
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
          print "   |!| please click the start of a sky region"
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
          print "Are you happy with the sky subtraction apertures? (default = no)"
          answer = raw_input("  (y)es, (n)o, (r)edo\n")
          if "y" in answer:
            finished = True
          elif "r" in answer:
            mask *= 0
          else:
            finished = False
        first = False
      self.images['skyMask'] = mask
      np.save(filename, self.images['skyMask'])
      print " Saved a sky mask to {0}".format(filename)

  def createWavelengthCal(self, visualize=True):
    '''Populate the wavelength calibration for this aperture.'''
    print "     populating wavelength calibration"
    filename = self.directory + 'waveCal_{0}.npy'.format(self.name)
    try:
      # load the the thing!
      self.waveCalCoef = np.load(filename)
    except:
      print "       estimating wavelength calibration from extracted arc spectra"
      # extract the arc lamp spectra
      self.arcs = {}
      for element in ['He', 'Ne','Ar']:
        self.arcs[element] = self.extract(self.images[element], arc=True)

      # load wavelength identifications
      wavelength_ids = astropy.io.ascii.read(self.obs.wavelengthFile)

      # cross correlate my arc spectra with reference
      def findRoughShift(wavelength_ids, visualize=False):
        print "       cross correlating arcs with known wavelengths"
        if visualize:
          fi, ax = plt.subplots(3,2, figsize=(5,5))

        count = 0
        for element in ['He', 'Ne','Ar']:
          flux = self.arcs[element]['raw_counts']
          x = self.arcs[element]['w']

          xPeak, yPeak = zachopy.oned.peaks(self.waxis, flux)

          this = []
          for i in range(len(wavelength_ids)):
            if element in wavelength_ids['name'][i]:
              print wavelength_ids[i]
              this.append(i)

          myPeaks, theirPeaks = np.zeros(len(x)), np.zeros(len(x))
          for i in range(len(this)):
            center = wavelength_ids['pixels'][this[i]]/2
            blob = 2
            theirPeaks += np.exp(-0.5*((x-center)/blob)**2)
          for i in range(len(xPeak)):
            # MY COORDINATES ARE REVERSED RELATIVE TO IRAFS!
            center = xPeak[i]
            blob = 2
            myPeaks += np.exp(-0.5*((x-center)/blob)**2)*np.log(yPeak[i])

          if visualize:
            ax[count,0].plot(x, theirPeaks/theirPeaks.sum())
            ax[count,0].plot(x, myPeaks/myPeaks.sum())
            ax[count,1].plot(np.correlate(myPeaks, theirPeaks, 'full'))

          if count == 0:
            corre = np.correlate(myPeaks, theirPeaks, 'full')
          else:
            corre = corre*np.correlate(myPeaks, theirPeaks, 'full')
          count += 1

        if visualize:
          for count in [0,1,2]:
            ax[count,1].plot(corre, color='black')
        peakoffset = np.where(corre == corre.max())[0][0] - len(x)
        # to convert: len(x) - xPeak = x + peakoffset

        return peakoffset

      # perform a first rough alignment, to make pixels
      peakoffset = findRoughShift(wavelength_ids)

      # empty lists, which we'll fill with pixels matched to wavelengths
      pixel, wavelength = [],[]

      # treat the arc lamps separately
      count =0
      for element in ['He', 'Ne','Ar']:

        # the pixel spectrum extracted from this arc lamp
        flux = self.arcs[element]['raw_counts']
        x = self.arcs[element]['w']

        xPeak, yPeak = zachopy.oned.peaks(self.waxis, flux)

        this = []
        for i in range(len(wavelength_ids)):
          if element in wavelength_ids['name'][i]:
            #print wavelength_ids[i]
            this.append(i)

        iMine, iTheirs =[],[]
        theirPeaks =  wavelength_ids['pixels'][this]/2
        theirWavelengths = wavelength_ids['wave1'][this]
        theirPeaksOnMyPixels = (theirPeaks +  peakoffset)
        myPeaksOnMyPixels = xPeak
        for i in range(len(theirPeaksOnMyPixels)):
          distance = xPeak - theirPeaksOnMyPixels[i]
          closestdistance = (np.abs(distance)).min()

          #print theirPeaksOnMyPixels[i], closestdistance, input['name'][this[i]]
          if closestdistance < 10:
            iMine.append(np.where(np.abs(distance) == closestdistance)[0][0])
            iTheirs.append(i)

        initialCoeff = np.polynomial.polynomial.polyfit(xPeak[iMine], theirWavelengths[iTheirs], 2)
        initialFit = np.polynomial.polynomial.Polynomial(initialCoeff)
        pixel.extend(xPeak[iMine])
        wavelength.extend(theirWavelengths[iTheirs])
        count += 1

      fit = np.polynomial.polynomial.Polynomial(np.polynomial.polynomial.polyfit(pixel, wavelength, 5) )
      residual = wavelength - fit(pixel)
      good = np.abs(residual) < np.median(np.abs(residual))*5
      bad = ~good
      self.waveCalCoef  = np.polynomial.polynomial.polyfit(pixel, wavelength, 5, w=good)
      fit = np.polynomial.polynomial.Polynomial(self.waveCalCoef )
      pixel, wavelength = np.array(pixel), np.array(wavelength)

      if visualize:
        plt.figure('wavelength calibration')
        ax = plt.subplot(211)
        ax.set_title('Wavelength Calib. for Aperture (%0.1f,%0.1f)' % (self.x, self.y))
        ax.scatter(pixel, wavelength, marker='o')
        ax.plot(pixel, fit(pixel))
        ax.set_ylabel('Wavelength (angstroms)')
        ax.set_ylabel('Residuals')
        ax = plt.subplot(212)
        ax.scatter(pixel[good], wavelength[good] - fit(pixel[good]), marker='o')
        ax.scatter(pixel[bad], wavelength[bad] - fit(pixel[bad]), marker='x', color='blue')
        ax.set_xlabel('Pixel # (by python rules)')
      np.save(filename, self.waveCalCoef)
      print "      saved wavelength calibration to ", filename
    self.wavelengthCalibrate = np.polynomial.polynomial.Polynomial(self.waveCalCoef )

  def ones(self):
    '''Create a blank array of ones to fill this aperture.'''
    return np.ones_like(self.images['Science'])

  def extract(self, raw, subtractsky=True, arc=False, cosmics=False, n=0):
    '''Extract the spectrum from this aperture.'''
    intermediates = {}
    intermediates['original'] = raw
    extracted = {}
    extracted['w'] = self.waxis

    # remove cosmic rays
    if cosmics and not arc:
      image = self.removeCosmics(raw)
    else:
      image = raw
    intermediates['cosmics'] = raw - image

    # define an image marking the distance from the ridge of the trace
    distanceFromTrace = np.abs(self.s - self.traceCenter(self.w))

    # define the extraction aperture (Gaussian profile for wavelength extraction, boxcar for stellar flux)
    if arc:
      intermediates['extractMask'] = self.images['RoughLSF']
      subtractsky = False
    else:
      intermediates['extractMask'] = self.ones()

    # replaced self.obs.nFWHM*width with "obs.extractionWidth"
    intermediates['extractMask'][distanceFromTrace > self.obs.extractionWidth] = 0

    # define a mask for the sky estimate area
    intermediates['skyMask'] = self.images['skyMask']
    #intermediates['skyMask'] = self.ones()
    #intermediates['skyMask'][distanceFromTrace < self.obs.extractionWidth + self.obs.skyGap] = 0
    #intermediates['skyMask'][distanceFromTrace > (self.obs.extractionWidth + self.obs.skyGap + self.obs.skyWidth)] = 0

    if subtractsky:
      # estimate a 1D sky spectrum by summing over the masked region
      #skyperpixel = (image*intermediates['skyMask']/self.images['NormalizedFlat']).sum(self.sindex)/(intermediates['skyMask']/self.images['NormalizedFlat']).sum(self.sindex)
      ma = np.ma.MaskedArray(image*intermediates['skyMask']/self.images['NormalizedFlat'], intermediates['skyMask']==0)
      skyperpixel = np.ma.median(ma, self.sindex).data

      extracted['sky'] = skyperpixel*intermediates['extractMask'].sum(self.sindex)
      intermediates['sky'] = np.ones_like(image)*skyperpixel.reshape((self.waxis.shape[0],1))

    # wavelength calibrate the spectrum, if you can
    try:
      extracted['wavelength'] = self.wavelengthCalibrate(self.waxis)
    except:
      extracted['wavelength'] = None

    # if sky subtraction is turned on, subtract off a sky spectrum
    if subtractsky:
      extracted['raw_counts'] = (intermediates['extractMask']*image/self.images['NormalizedFlat']).sum(self.sindex) - extracted['sky']
      intermediates['subtracted'] = np.maximum(image - intermediates['sky'], 0)
      extracted['centroid'] = np.nansum(self.s*intermediates['subtracted']*intermediates['extractMask'],self.sindex)/np.nansum(intermediates['subtracted']*intermediates['extractMask'], self.sindex)
      extracted['width'] = np.sqrt(np.nansum(self.s**2*intermediates['subtracted']*intermediates['extractMask'], self.sindex)/np.nansum(intermediates['subtracted']*intermediates['extractMask'], self.sindex) - extracted['centroid']**2)
      extracted['peak'] = np.max(intermediates['extractMask']*image, self.sindex)
      #assert(np.isfinite(extracted['centroid']).all())

    else:
      extracted['raw_counts'] = np.nansum(intermediates['extractMask']*image/self.images['NormalizedFlat'], self.sindex)

    #self.displayStamps(intermediates)
    #if subtractsky:
      #for k in ['subtracted']:# intermediates.keys():
      #	writeFitsData(intermediates[k], self.directory + k + '{0:04}.fits'.format(n))

    if self.visualize:
      self.display.rgb(intermediates['original'], intermediates['extractMask'], intermediates['skyMask'])
      a = raw_input("OK?")
      if 'y' in a:
        self.visualize=False
    return extracted

  def plot(self, extracted, coordinate='wavelength', sharex=True, filename=None):
    '''Plot extracted spectrum.'''

    # switch to wavelength coordinates, if necessary (and possible)
    x = extracted['w']
    if coordinate=='wavelength':
      x = extracted['wavelength']

    # select the large figure, and clear it
    keys = extracted.keys()
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
      self.eax[i].plot(x,extracted[keys[i]])
      self.eax[i].set_ylabel(keys[i])
    self.eax[-1].set_xlim(x.min(), x.max())
    self.eax[-1].set_xlabel(coordinate)
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

  def extractAll(self, remake=False):
    '''Loop through exposures, extracting all the spectra.'''

    for n in self.obs.nScience:
      extractedFilename = self.directory + 'extracted{0:04}.npy'.format(n)
      if os.path.exists(extractedFilename) and remake == False:
        print "       {0} already extracted".format(extractedFilename)
      else:
        try:
          stamp = readFitsData(stampFilename(n,self), verbose=True)
        except:
          self.mask.createStamps(n)
          stamp = readFitsData(stampFilename(n,self), verbose=True)
        extracted = self.extract(stamp, n=n)
        np.save(extractedFilename, extracted)
        if self.visualize:
          self.plot(extracted)
