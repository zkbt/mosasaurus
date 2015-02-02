from imports import *

class Cube():
  '''Cube object stores a -wavelength-flux datacube.'''
  def __init__(self, obs, remake=False):
    '''Initialize a data cube and populate it with data.'''
    self.obs = obs
    self.display = self.obs.display
    self.tempfilename = self.obs.extractionDirectory + 'tempSpectralCube.npy'
    self.filename = self.obs.extractionDirectory + 'spectralCube.npy'
    self.keys = ['sky',  'centroid', 'width', 'peak', 'raw_counts']
    self.loadSpectra(remake=remake)
    self.sindex = 0
    self.tindex = 1
    self.windex = 2
    self.goodComps = self.obs.goodComps

  def loadSpectra(self, remake=False, visualize=False):
    """Opens all the spectra in a working directory, lumps them into a cube, and returns it:
      cube, wave = LDSS3.loadSpectra(remake=True)
      cube = an array containing extracted spectra of shape (nstars, nexposures, nwavelengths)
      wave = an array containing wavelength of shape (nwavelengths)

      by default, will search for a stored "cube.npy" file in working directory to load;
        if change have been made, run with remake=True"""
    print "   Loading the spectral cube."
    try:
      assert(remake == False)
      self.commonwavelength,self.dnpixelsdw,self.cubes = np.load(self.filename)
    except:


      # define the directories that contain extracted stellar spectra
      starDirectories = glob.glob(self.obs.extractionDirectory + 'aperture_*')

      # set up a fine, common wavelength grid onto which everything will be interpolated
      commonwavelength = np.arange(4000, 10500)

      # set up a plot window to show how the interpolation is going
      self.display.medium()

      # define an empty cube that we're going to populate with spectra
      self.cubes = {}
      if visualize:
        self.ax_cube = {}
        gs = plt.matplotlib.gridspec.GridSpec(len(keys),1,hspace=0,wspace=0)





      numberofstars = len(starDirectories)
      numberofspectra = len(self.obs.nScience)

      # loop over all the stars
      for star in range(numberofstars):



        # loop over the spectra
        for spectrum in range(numberofspectra):

          # find the available spectra
          spectrumFile = starDirectories[star] + '/extracted{0:04.0f}.npy'.format(self.obs.nScience[spectrum])
          print "   ", spectrumFile

          # load the extracted spectrum
          extracted = np.load(spectrumFile)[()]
          wavelength = extracted['wavelength']

          try:
            self.commonwavelength
          except:
            # calculate the number of pixels that go into each wavelength bin
            dw_original = wavelength[1:] - wavelength[0:-1]
            dw_new = np.ones_like(commonwavelength)
            interpolation = interp.interp1d(wavelength[0:-1], dw_original, bounds_error=False)
            dnpixelsdw = dw_new/interpolation(commonwavelength)
            self.commonwavelength = commonwavelength
            self.dnpixelsdw = dnpixelsdw

          # loop over the measurements
          count = 0
          for key in self.keys:
            if star == 0 and spectrum == 0:

              # create an empty cube
              self.cubes[key] = np.zeros((numberofstars, numberofspectra, len(commonwavelength)))

              # set up the plots
              if visualize:
                if key != 'raw_counts':
                  sharex = self.ax_cube['raw_counts']
                else:
                  sharex = None
                self.ax_cube[key] = plt.subplot(gs[count], sharex=sharex)
                count += 1
            else:
              # clear the plots
              if visualize:
                self.ax_cube[key].cla()


            # supersample onto the grid
            self.cubes[key][star,spectrum,:] = zachopy.oned.supersample(wavelength, extracted[key], commonwavelength, visualize=False)

            # plot demonstration
            if visualize and True:
              self.ax_cube[key].set_ylabel(key)
              self.ax_cube[key].plot(wavelength, extracted[key], color='black', alpha=0.5)
              self.ax_cube[key].plot(commonwavelength, self.cubes[key][star,spectrum,:], color='red', alpha=0.5)

              if key != keys[-1]:
                plt.setp(self.ax_cube[key].get_xticklabels(), visible=False)

          if visualize:
            plt.draw()
          #a = raw_input("supersampling okay?")

          # debug!
          if spectrum % 100 == 0:
            print "saved (temporarily) to ", self.tempfilename
            np.save(self.tempfilename,(self.commonwavelength,self.dnpixelsdw,self.cubes))

      np.save(self.filename,(self.commonwavelength,self.dnpixelsdw,self.cubes))
    print "   Done loading spectral cube."

  def convolveCube(self,width=0.5):
    '''Take the cube of spectra, and convolve the spectra by a Gaussian.'''
    c = self.flux
    w = self.commonwavelength
    '''Take a spectral cube, convolve it with a Gaussian.'''
    star = self.obs.target
    start = zachopy.oned.subtractContinuum(c[star,0,:])
    plt.ion()
    x = np.arange(-width*10, width*10)
    gauss = lambda x: np.exp(-0.5*x**2/width**2)
    i = 0
    new_cube = np.copy(c)
    for star in range(len(c[:,0,0])):
      plt.figure()
      plt.plot(c[star,i,:])
      plt.plot(np.convolve(c[star,i,:], gauss(x)/np.sum(gauss(x)),'same'))
      for i in range(len(c[0,:,0])):
        new_cube[star,i,:] = np.convolve(c[star,i,:], gauss(x),'same')
    self.flux = new_cube

  def shiftCube(self,plot=True):
    '''Shift all the spectra for a star to line up in wavelength.'''
    star = self.obs.target
    c = self.flux
    '''Take a spectral cube, shift all the spectra to line up with one.'''
    plt.figure()
    plt.ion()
    left = np.argmin(np.abs(self.commonwavelength - 8350))
    right = np.argmin(np.abs(self.commonwavelength - 8800))
    width = 5
    triplet = [8498.,8542.,8662.]
    # shift all O lines to one star!
    wave = self.commonwavelength[left:right]
    start = np.ones_like(wave)*np.std(zachopy.oned.subtractContinuum(c[self.goodComps[0],c.shape[1]/2,left:right]))*3
    for t in triplet:
      start *= 1 - np.exp(-(wave - t)**2/2/width**2)
    start -= np.median(start)
    new_cube = np.copy(c)
    star_shifts = []

    for star in range(len(c[:,0,0])):
      shifts = []
      for i in range(len(c[0,:,0])):
        if i > -1:
          this = zachopy.oned.subtractContinuum(c[star,i,left:right])
          xc = np.correlate(this,start,'same')
          x = np.arange(len(xc))
          coeff = np.polynomial.polynomial.polyfit(x[xc.argmax()-5:xc.argmax()+5], xc[xc.argmax()-5:xc.argmax()+5], 2)
          fit = np.polynomial.polynomial.Polynomial(coeff)
          der = fit.deriv()
          peak = der.roots()
          #plt.plot(x, xc, alpha=0.3)
          #plt.scatter(peak, i)

          offset = peak - len(start)/2 + 1

          pixels = np.arange(len(c[star,i,:]))
          interpolation = interp.interp1d(pixels, c[star,i,:], kind='linear', bounds_error=False, fill_value=0.0)
          new_cube[star,i,:] = interpolation(pixels + offset)
          if plot:
            if i == 0:
              fi,ax = plt.subplots(2,1)
            else:
              for a in ax:
                a.cla()
            ax[0].plot(x, xc, alpha=0.3)
            #ax[0].plot(x,fit(x))
            ax[0].scatter(peak, fit(peak))
            ax[0].set_xlim(0,len(x))
            ax[1].plot(wave, start, color='black')
            ax[1].plot(wave, start, color='black', alpha=0.2, linewidth=5)
            ax[1].plot(wave, this, color='gray', alpha=0.2, linewidth=5)
            new = zachopy.oned.subtractContinuum(new_cube[star,i,left:right])
            ax[1].plot(wave, new, color='green', alpha=0.9, linewidth=2)
            ax[1].set_xlim(wave.min(), wave.max())
            ax[1].set_autoscaley_on
            ax[0].set_title('{0}/{1} stars; {2}/{3} spectra'.format(star+1, len(c[:,0,0]), i+1, len(c[0,:,0])))
            ax[0].axvline(len(x)/2.0)
            plt.draw()

          shifts.append(offset)
          print "shift = {4} for {0}/{1} stars; {2}/{3} spectra".format(star+1, len(c[:,0,0]), i+1, len(c[0,:,0]), offset)
      star_shifts.append(shifts)
    self.flux = new_cube
    self.shifts = star_shifts

  def createBins(self, binsize=250, remake=False):
    '''Take a spectral cube, and creates big wavelength bins.'''
    binned_filename = self.obs.extractionDirectory + "binnedcube{0:04}.npy".format(binsize)
    try:
       (self.binned_cubes, self.bin_centers, self.bin_ok) = np.load(binned_filename)
       assert(remake==False)
       print "   A binned cube for {0}A bins already exists. Loaded it!".format(binsize)
    except:
      print "   Binning the spectral cube."
      #self.shiftCube(plot=True)
      #self.show()
      wavelength = self.commonwavelength
      dnpixelsdw = self.dnpixelsdw

      # the shape of the cube
      nStars, nTimes, nWaves = self.cubes['raw_counts'].shape
      plt.ion()

      self.binned_cubes = {}


      # define bins
      bin_starts = np.arange(wavelength.min(), wavelength.max() - binsize/2.0, binsize)
      bin_ends = bin_starts + binsize
      bin_centers = (bin_starts + bin_ends)/2.0
      nBins = bin_centers.size
      bin_ok = np.zeros((nStars, nBins))
      satlimit = 150000
      faintlimit = 2500

      # loop through stars, determining which bins are okay for each
      for star in range(nStars):
        for bin in range(nBins):
          mask = ( wavelength >= bin_starts[bin] ) * (wavelength <= bin_ends[bin])
          notsaturated = np.median(np.max(self.cubes['peak'][star,:,mask], 0) < satlimit)
          brightenough = np.median(np.max(self.cubes['raw_counts'][star,:,mask], 0) > faintlimit)
          bin_ok[star, bin] = notsaturated*brightenough
          print '-------------'
          print '   star {0} is {1}'.format(star, bin_ok[star,:])
        if star == 1:
          assert(bin_ok[star,:].any())

      # create a binned cube
      for key in self.keys:
        if key == 'flux':
          newkey = 'raw_counts'
        else:
          newkey = key
        self.binned_cubes[newkey] = np.ones((nStars, nTimes, nBins)).astype(np.float)

      # loop over the stars
      for star in range(nStars):
        # loop over the times
        for time in range(nTimes):
          print "         star {0:10}, timepoint {1:10}".format(star, time)
          # loop over the bins
          for bin in range(nBins):

            # loop over keys
            for key in self.keys:
              if key == 'flux':
                newkey = 'raw_counts'
              else:
                newkey = key
              wholespectrum = self.cubes[key][star,time,:]
              # pick out the wavelengths that are on this bin
              mask = (wavelength > bin_starts[bin])*(wavelength < bin_ends[bin])

              # integrate spectrum over bin (dndw factor accounts for the fractional pixels represented by resampled cube spectrum
              self.binned_cubes[newkey][star, time, bin] = scipy.integrate.trapz(wholespectrum[mask]*dnpixelsdw[mask], wavelength[mask])

      self.bin_centers = bin_centers
      self.bin_ok = bin_ok
      np.save(binned_filename, (self.binned_cubes, self.bin_centers, self.bin_ok))
      print "    Saved binned cube to {0}".format(binned_filename)

  def imageBins(self, binned_cube, title=' ', vmin=None, vmax=None):
    '''Make an image of the input cube.'''
    print "   Trying to image bins for " + title
    bin_centers = self.bin_centers
    bin_ok = self.bin_ok

    '''Take binned lightcurves, image them. Returns the bins with median correction subtracted.'''
    plt.ion()
    if vmin == None:
      vmin = 0.98
    if vmax == None:
      vmax = 1.01
    nStars, nTimes, nWaves = binned_cube.shape
    fi, ax = plt.subplots(nStars, 1, figsize=(10,10), sharex=True, sharey=True)
    fi.suptitle(title)
    binsize = bin_centers[1] - bin_centers[0]
    self.bin_starts = bin_centers - binsize/2.0
    self.bin_ends = bin_centers + binsize/2.0
    nBins = bin_centers.size

    corrected = np.zeros_like(binned_cube)
    for star in range(nStars):
      ax[star].imshow(binned_cube[star,:,:], vmin=vmin, vmax=vmax, extent=[self.bin_starts.min(), self.bin_ends.max(), 0, nTimes], aspect='auto', cmap='gray', interpolation='nearest')
      for bin in range(nBins):
        if bin_ok[star,bin] == 0:
          ax[star].axvspan(self.bin_starts[bin], self.bin_ends[bin], alpha=0.7, color='white', edgecolor=None)
    a = raw_input('imaging ' + title + '?')

  def correctBins(self):
    '''Use comparison stars to correct for atmospheric losses, create a self.binned_corrected.'''


    wavelength = self.bin_centers
    ccdn, times, airmass = np.loadtxt(self.obs.workingDirectory + 'headerInfo.txt', unpack=True)
    vmin = 0.98
    vmax = 1.02
    nStars, nTimes, nWaves = self.binned_cubes['raw_counts'].shape

    # define the bins
    bin_centers = self.bin_centers
    binsize = bin_centers[1] - bin_centers[0]
    bin_starts = bin_centers - binsize/2.0
    bin_ends = bin_centers + binsize/2.0
    nBins = bin_centers.size

    correction = np.ones((nTimes, nWaves))
    corrected = np.ones_like(self.binned_cubes['raw_counts'])

    for star in self.goodComps:
      for time in range(nTimes):
        mask = np.bool8(self.bin_ok[star,:])
        correction[time,mask] += self.binned_cubes['raw_counts'][star,time,mask]

    # normalize the correction spectrum to be close to one
    mediancompositespectrum = np.median(correction, 0)
    for time in range(nTimes):
      correction[time,:] /= mediancompositespectrum
    self.binned_cubes['corrected'] = self.binned_cubes['raw_counts']/correction
    self.binned_cubes['error'] = np.sqrt(self.binned_cubes['raw_counts'] + self.binned_cubes['sky'])

  def fitOot(self):
    '''Subtract off a polynomial fit to out of transit continuum.'''
    cube = self.binned_corrected
    nStars, nTimes, nWaves = cube.shape
    mask = np.zeros(nTimes)
    #mask[108:185] = True
    #mask[55:140] = True
    mask[250:900] = True
    t = np.arange(nTimes)
    plt.ion()
    plt.figure()
    clean_cube = np.zeros_like(cube)
    for star in range(nStars):
      for wave in range(nWaves):
        timeseries = cube[star,:,wave]
        fit = np.polynomial.polynomial.Polynomial(np.polynomial.polynomial.polyfit(t[mask==False], timeseries[mask==False], 2))
        clean_cube[star,:,wave] = cube[star,:,wave]/fit(t)
    self.binned_cleaned = clean_cube

  def plotBins(self):
    cleaned = self.binned_cleaned
    bin_centers = self.bin_centers
    bin_ok = self.bin_ok

    '''Take some binned lightcurves, plot them.'''
    star = self.obs.target
    nStars, nTimes, nWaves = cleaned.shape
    fi, ax = plt.subplots(int(np.sqrt(nWaves))+1, int(np.sqrt(nWaves))+1, sharex=True, sharey=True, figsize=(10,10))
    flat = ax.flatten()
    for wave in range(nWaves):
      if bin_ok[star,wave]:
        alpha = 1.0
      else:
        alpha = 0.3
      flat[wave].plot(cleaned[star,:,wave], alpha=alpha)
      flat[wave].set_ylim(0.97, 1.01)
      flat[wave].set_xlim(0, nTimes)
      flat[wave].set_title('%d' % bin_centers[wave])

  def wasp94(self, binsize=200):
    '''Kludge for plotting some WASP-94 rough telescope data.'''
    star = self.obs.target

  '''def imagingtemp(self):

    # show an image of the raw light curves
    binned_cube = self.binned_cubes['raw_counts']
    normalized = np.ones_like(binned_cube)
    for star in np.arange(binned_cube.shape[0]):
      medianspectrum = np.median(binned_cube[star,:,:],0)*np.ones_like(binned_cube[star,:,:])
      normalized[star, :, :] = binned_cube[star, :, :]/medianspectrum
    self.imageBins(normalized, title='Uncorrected')


    binned_cube = self.binned_corrected
    normalized = np.ones_like(binned_cube)
    for star in np.arange(binned_cube.shape[0]):
      medianspectrum = np.median(binned_cube[star,:,:],0)*np.ones_like(binned_cube[star,:,:])
      normalized[star, :, :] = binned_cube[star, :, :]/medianspectrum
    self.imageBins(normalized, title='Corrected')

    # for each bin, correct lightcurves by fitting polynomial to oot
    #self.fitOot()
    #self.imageBins(self.binned_cleaned, title='Cleaned')

    '''
  def makeMeanSpectrum(self):
    self.loadSpectra()
    wavelength = self.commonwavelength
    spectrum = np.median(self.cubes['raw_counts'][self.obs.target,:,:],1).flatten()
    assert(len(spectrum) == len(wavelength))
    fi, ax = plt.subplots(1)
    unit = 10.0
    ax.plot(wavelength/unit, spectrum*unit)
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Flux (photons/nm/exposure)')
    print "   Saving median spectrum to " + self.obs.extractionDirectory + 'medianSpectrum.npy'
    np.save(self.obs.extractionDirectory + 'medianSpectrum.npy', (wavelength, spectrum))

  def makeLCs(self,binsize=250):
    '''Wrapper to go from extracted spectra to binned, multiwavelength lightcurves.'''

    # pick the target star
    star = self.obs.target

    # bin the cube into manageable wavelength bins
    self.createBins(binsize=binsize)


    # use comparison star(s) to divide out flux losses
    self.correctBins()




    # setup
    nStars, nTimes, nWaves = self.binned_cubes['corrected'].shape
    bw = self.bin_centers
    binsize = bw[1] - bw[0]
    bin_starts = bw - binsize/2.0
    bin_ends = bw + binsize/2.0



    lcDirectory = self.obs.extractionDirectory + 'lc_binby' + ('%d' % binsize) + '/'
    if not os.path.exists(lcDirectory):
      os.makedirs(lcDirectory)

    try:
      ccdn, times, airmass = np.loadtxt(self.obs.workingDirectory + 'headerInfo.txt', unpack=True)
      assert(len(times) == self.binned_cubes['raw_counts'].shape[1])
    except:
      self.obs.loadHeaders()
      ccdn, times, airmass = np.loadtxt(self.obs.workingDirectory + 'headerInfo.txt', unpack=True)

    try:
      idl = pidly.IDL('/Applications/exelis/idl/bin/idl')
      bjd = idl.func('utc2bjd', times + 2400000.5, self.obs.ra, self.obs.dec)
    except:
      print "====================================="
      print "====================================="
      print "========   NO BJD CORRECTION! ======="
      print "====================================="
      print "====================================="
      print "====================================="
      bjd = times + 2400000.5

    self.lcs = []
    for wave in range(nWaves):
      target = self.obs.target
      if self.bin_ok[target,wave]:
        if np.isfinite(self.binned_cubes['corrected'][target,:,wave]).any():
          lc = LC(self.obs, bin_starts[wave], bin_ends[wave])
          flux = self.binned_cubes['corrected'][target,:,wave].flatten()/np.median(self.binned_cubes['raw_counts'][target,:,wave].flatten())
          error = self.binned_cubes['error'][target,:,wave].flatten()/self.binned_cubes['raw_counts'][target,:,wave].flatten()
          dict = {}
          for key in self.binned_cubes.keys():
            if key != 'corrected' and key != 'error':
              dict[key+'_target'] = self.binned_cubes[key][target,:,wave].flatten()
              for comparison in self.obs.goodComps:
                dict[key+'_comparison{0:02.0f}'.format(comparison)] = self.binned_cubes[key][comparison,:,wave].flatten()
          dict['airmass'] = airmass
          lc.populate(bjd, flux, error, **dict)

          #lc.plot()
          lc.save()
          self.lcs.append(lc)

  def loadLCs(self, binsize=250):
    lcDirectory = self.obs.extractionDirectory + 'lc_binby' + ('%d' % binsize) + '/'
    g = glob.glob(lcDirectory + 'lc_*.npy')
    wavelengths = []
    lcs = []
    for file in g:
      lc = LC(self.obs, filename=file)
      if lc.lc is not None:
        wavelengths.append(lc.wavelength)
        lcs.append(lc)
    self.lcs = np.array(lcs)
    self.lcs = self.lcs[np.argsort(wavelengths)]
    return self.lcs

  def show(self, binsize=500,vmin = 0.98,vmax = 1.02):

      cube = self.flux
      wavelength = self.commonwavelength
      plt.ion()
      ccdn, times, airmass = np.loadtxt(self.obs.workingDirectory + 'headerInfo.txt', unpack=True)
      ok = np.arange(len(cube[0,:,0]))
      ccdn = ccdn[ok]
      times = times[ok]
      airmass = airmass[ok]


      nStars, nTimes, nWaves = cube.shape
      fi, ax = plt.subplots(nStars, 2, figsize=(10,10), sharex=True, sharey=True)
      fi.subplots_adjust(hspace=0, wspace=0)
      dividedbymedianspectrum = cube/np.median(cube,self.tindex).reshape(cube.shape[0], 1, cube.shape[2])
      for i in range(nStars):
        ax[i,0].imshow(dividedbymedianspectrum[i,:,:], vmin=vmin, vmax=vmax, extent=[wavelength.min(), wavelength.max(), 0, nTimes], aspect='auto',  cmap='gray', interpolation='nearest')

      sum = np.sum(cube[self.goodComps,:,:], self.sindex)
      sum /= np.median(sum,0).reshape(1, sum.shape[1])
      dividedbycomps =  dividedbymedianspectrum/sum.reshape(1, cube.shape[1], cube.shape[2])

      #	ax[1,1].imshow(sum, vmin=vmin, vmax=vmax, extent=[wavelength.min(), wavelength.max(), 0, nTimes], aspect='auto')
      for i in range(nStars):
        ax[i,1].imshow(dividedbycomps[i,:,:], vmin=vmin, vmax=vmax, extent=[wavelength.min(), wavelength.max(), 0, nTimes], aspect='auto', cmap='gray', interpolation='nearest')

      bin_starts = np.arange(wavelength.min(), wavelength.max() - binsize/2.0, binsize)
      bin_ends = bin_starts + binsize
      bins = (bin_starts + bin_ends)/2.0
      nBins = bins.size

      fi, ax =  plt.subplots(len(bins) + 2, nStars, figsize=(15,10))
      fi.subplots_adjust(hspace=0, wspace=0)

      for i in range(nStars):
          ax[0,i].plot(wavelength, cube[i,:,:].sum(0))
          ax[0,i].set_title('%d' % i)
          ax[-1,i].plot(times, airmass)

      row = 0
      for b in bins:
        row += 1
        waveMin = b - binsize/2
        waveMax = b + binsize/2
        iWave = (wavelength >= waveMin) * (wavelength <= waveMax)
        #print b, ': ', wavelength[iWave]

        for i in range(nStars):
          ax[row,i].plot(times, np.nansum((cube[i,:,iWave].transpose()/sum[:,iWave]), 1))

      for axis in ax.reshape(-1):
        axis.xaxis.set_ticklabels([])
        axis.yaxis.set_ticklabels([])

  def divide(self):
    # calculate a median transit shape, and divide through
    self.target_raw = self.binned_cube[self.obs.target[0],:,:]
    self.target_corrected = self.binned_corrected[self.obs.target[0],:,:]
    transit = scipy.signal.medfilt(np.median(self.target_corrected, -1),3)

    self.target_median = transit.reshape(transit.shape[0], 1)*np.ones_like(self.target_corrected)
    self.target_divided = self.target_corrected/self.target_median


    plt.figure()
    plt.plot(transit)
    self.imageTarget()




  def imageTarget(self, title=None, vmin=None, vmax=None):
    '''Make an image of the input cube.'''
    if title is None:
      title = self.obs.name + ' | ' + self.obs.night
    print "   Trying to image bins for " + title
    bin_centers = self.bin_centers
    bin_ok = self.bin_ok
    target = self.target_raw
    '''Take binned lightcurves, image them. Returns the bins with median correction subtracted.'''
    plt.ion()
    if vmin == None:
      vmin = 0.985
    if vmax == None:
      vmax = 1.005
    nTimes, nWaves = target.shape
    fi, ax = plt.subplots(4,1, figsize=(4.5,12), sharex=True, sharey=True)
    plt.subplots_adjust()
    ax[0].set_title(title)
    binsize = bin_centers[1] - bin_centers[0]
    self.bin_starts = bin_centers - binsize/2.0
    self.bin_ends = bin_centers + binsize/2.0
    nBins = bin_centers.size

    normalized = np.ones_like(self.target_raw)
    targets = [self.target_raw, self.target_corrected, self.target_median, self.target_divided]
    names = ['Raw Detected Flux', 'Divided by Comparison', 'Median Transit', 'Divided by Transit']
    for i in range(len(targets)):
      target = targets[i]

      medianspectrum = np.median(target,0)*np.ones_like(target)
      normalized = target/medianspectrum
      ax[i].imshow(normalized, vmin=vmin, vmax=vmax, extent=[self.bin_starts.min(), self.bin_ends.max(), 0, nTimes], aspect='auto', cmap='gray', interpolation='nearest')
      ax[i].set_ylabel(names[i])
    ax[-1].set_xlabel('Wavelength (angstroms)')
    ax[-1].set_xlim(self.bin_starts.min(), self.bin_ends.max())
    for bin in range(nBins):
      if bin_ok[self.obs.target,bin] == 0:
        for a in ax:
          a.axvspan(self.bin_starts[bin], self.bin_ends[bin], alpha=0.7, color='white', edgecolor=None)
    plt.tight_layout(h_pad=0.0)
    plt.savefig(self.obs.workingDirectory + 'divided_{0}.pdf'.format(self.obs.night))
    a = raw_input('imaging  target')

class LC():
  def __init__(self, obs, left=None, right=None, filename=None):
    self.obs = obs
    if filename is not None:
      self.load(filename)
    else:
      self.left = left
      self.right = right
      self.setup()


  def setup(self):
    self.wavelength = (self.left + self.right)/2.0
    self.binsize = self.right - self.left
    self.filename = self.obs.extractionDirectory + 'lc_binby' + ('%d' % self.binsize) + '/lc_{0:05d}to{1:05d}.npy'.format(np.int(self.left), np.int(self.right))

  def populate(self, bjd, flux, error, **kwargs):
    # set up the column names for the light curve record array
    types = [('bjd', np.float), ('flux', np.float), ('error', np.float)]
    for key in kwargs.keys():
      types.append((key,np.float))

    # populate the columns with data
    self.lc = np.zeros(len(bjd), types)
    self.lc['bjd'] = bjd
    self.lc['flux'] = flux
    self.lc['error'] = error
    for key in kwargs.keys():
      self.lc[key] = kwargs[key]



  def save(self):
    np.save(self.filename, self.lc)

  def load(self, filename):
    self.left = np.float(filename.split('lc_')[-1].split('to')[-2])
    self.right = np.float(filename.split('.npy')[-2].split('to')[-1])
    self.setup()
    assert(self.filename == filename)
    self.lc = np.load(self.filename)


  def plot(self):
    try:
      for a in self.ax.values():
        a.cla()
    except:
      self.keystoplot = ['flux', 'raw_counts_target', 'sky_target', 'airmass', 'width_target', 'centroid_target', 'peak_target']
      gs = plt.matplotlib.gridspec.GridSpec(len(self.keystoplot), 1, height_ratios=[5,1,1,1,1,1,1], wspace=0, hspace=0)
      self.ax = {}
      for i in range(len(self.keystoplot)):
        key = self.keystoplot[i]
        try:
          sharex = self.ax[self.keystoplot[0]]
        except:
          sharex = None
        self.ax[key] = plt.subplot(gs[i], sharex=sharex)
    for key in self.keystoplot:
      self.ax[key].plot(self.lc['bjd'], self.lc[key], markersize=1, marker='o', alpha=0.3, color='black', linewidth=0)
      self.ax[key].set_xlim(self.lc['bjd'].min(), self.lc['bjd'].max())
      self.ax[key].set_ylabel(key)

    plt.draw()
    bla = raw_input(self.filename + '?')
