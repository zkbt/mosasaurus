from imports import *
from Observation import Observation
import astropy.table, astropy.time
import pidly

plt.ion()
class Cube(Talker):
  '''Cube object stores a -wavelength-flux datacube.'''
  def __init__(self, obs, remake=False, max=None, **kwargs):
    '''Initialize a data cube and populate it with data.'''
    Talker.__init__(self, line=200, **kwargs)

    if type(obs) == str:
        self.obs = Observation(obs, nods9=True)
    else:
        self.obs = obs

    self.tempfilename = self.obs.extractionDirectory + 'tempSpectralCube.npy'
    self.cubekeys = ['sky',  'centroid', 'width', 'peak', 'raw_counts']

    self.savable = ['cubes', 'squares', 'temporal', 'spectral']
    self.sindex = 0
    self.tindex = 1
    self.windex = 2
    self.goodComps = self.obs.goodComps


  def populate(self, remake=False, max=None):
    try:
        self.cubes
    except:
        try:
            assert(remake==False)
            self.load()
        except:
            self.loadSpectra(remake=remake, max=max)

  @property
  def display(self):
    try:
        return self._display
    except:
        self._display = zachopy.display.ds9('cube')
        return self._display

  @property
  def filename(self):
    return self.obs.extractionDirectory + 'spectralCube_{0}stars_{1}spectra.npy'.format(self.numberofstars, self.numberoftimes)

  def loadSpectra(self, remake=False, visualize=True, max=None):
    """Opens all the spectra in a working directory, lumps them into a cube, and returns it:
      cube, wave = LDSS3.loadSpectra(remake=True)
      cube = an array containing extracted spectra of shape (nstars, nexposures, nwavelengths)
      wave = an array containing wavelength of shape (nwavelengths)

      by default, will search for a stored "cube.npy" file in working directory to load;
        if change have been made, run with remake=True"""

    # 3d, stars x time x wavelength
    self.cubes = {}
    # 2d, stars x time
    self.squares = {}
    # 1d, time
    self.temporal = {}
    # 1d, wavelength
    self.spectral = {}

    # update
    self.speak("Loading the spectral cube.")
    # define the directories that contain extracted stellar spectra
    starDirectories = glob.glob(self.obs.extractionDirectory + 'aperture_*')
    self.numberofstars = 2#len(starDirectories)
    self.numberoftimes = len(self.obs.nScience)
    if max is not None:
        self.numberoftimes = max
    truncate = False

    try:
        assert(remake==False)
        self.load()
        return
    except:
        self.speak('Could not load pre-saved cube. Creating one now!')
    # load the headers
    self.headers = astropy.table.Table(np.load(self.obs.workingDirectory + 'headers.npy')[()])

    # loop over the spectra
    for spectrum in range(self.numberoftimes):
        # loop over all the stars
        for star in range(self.numberofstars):

          # ccd number for this image
          ccdn =  self.obs.nScience[spectrum]

          # find the available spectrum
          spectrumFile = starDirectories[star] + '/supersampled{0:04.0f}.npy'.format(ccdn)
          extractedFile = starDirectories[star] + '/extracted{0:04.0f}.npy'.format(ccdn)

          self.speak('trying to load {0}'.format(spectrumFile))

          # load the extracted spectrum (or truncate the cubes at this point)
          try:
              supersampled = np.load(spectrumFile)[()]
              self.speak('loaded {0}'.format(spectrumFile))
              extracted = np.load(extractedFile)[()]
              self.speak('loaded {0}'.format(extractedFile))
              #print supersampled.keys()
              #print extracted.keys()
          except:
              # if we've run out of spectra to load, then truncate
              truncate = True
              breakf

          # loop over the measurement types and populate the cubes
          for key in self.cubekeys:
            try:
                self.cubes[key]
            except:
                self.cubes[key] = np.zeros((self.numberofstars, self.numberoftimes, len(supersampled['wavelength'])))
            self.speak('updating cubes[{0}][{1},{2},:]'.format(key,star, spectrum))
            # supersample onto the grid
            self.cubes[key][star,spectrum,:] = supersampled[key]

          # pull out data from the (unsupersampled) spectra to populate a square with dimensions self.numberofstars x self.numberoftimes
          ok = (extracted['wavelength'] > np.min(supersampled['wavelength']))*(extracted['wavelength'] < np.max(supersampled['wavelength']))
          for key in ['sky', 'width', 'centroid', 'cosmicdiagnostic']:
              try:
                  self.squares[key]
              except:
                  self.squares[key] = np.zeros((self.numberofstars, self.numberoftimes))
              self.squares[key][star,spectrum] = np.median(extracted[key])
              self.speak('updating squares[{0}][{1},{2}]'.format(key,star,spectrum))

        # if we've run out of spectra, then break out of the loop (with truncated cubes)
        if truncate:
            break


        self.speak('{0}/{1} spectra loaded into cube'.format(spectrum, self.numberoftimes))
        # if the spectra for all stars were successfully loaded, then
        try:
            self.temporal['n']
        except:
            self.temporal['n'] = []
        self.temporal['n'].append(ccdn)

    # make sure everything is truncated properly
    if truncate:
        self.speak("couldn't find all requested spectra, so truncated cube at a length of {0}".format(spectrum))
        for key in self.cubes.keys():
            self.cubes[key] = self.cubes[key][:,0:spectrum,:]
        for key in self.squares.keys():
            self.squares[key] = self.squares[key][:,0:spectrum]

    satlimit = 150000
    faintlimit = 2500
    self.cubes['ok'] = (self.cubes['peak'] < satlimit)*(self.cubes['raw_counts'] > faintlimit)
    self.cubes = astropy.table.Table(self.cubes)
    self.squares = astropy.table.Table(self.squares)
    self.temporal = astropy.table.join(self.headers, astropy.table.Table(self.temporal))
    self.temporal['cosmics'] = self.squares['cosmicdiagnostic'].sum(0)
    self.temporal['width'] = self.squares['width'][self.obs.target[0],:]
    self.temporal['centroid'] = self.squares['centroid'][self.obs.target[0],:]
    self.temporal['sky'] = self.squares['sky'][self.obs.target[0],:]
    self.temporal['lc'] = self.cubes['raw_counts'][0,:,:].sum(-1)/self.cubes['raw_counts'][1,:,:].sum(-1)

    self.temporal['ok'] = self.temporal['cosmics'] < 200

    self.cubes['ok'] *= self.temporal['ok'].reshape(1,self.numberoftimes,1)

    savable = ['cubes', 'squares', 'temporal', 'spectral']
    assert(len(self.cubes[0][0]) == len(self.squares[0][0]))
    assert(len(self.cubes[0][0]) == len(self.temporal))



    # define some useful arrays
    self.spectral['wavelength'] = supersampled['wavelength']
    self.spectral['dnpixelsdw'] = supersampled['dnpixelsdw']
    self.numberofwavelengths = len(self.spectral['wavelength'])

    self.shiftCube()
    self.speak("Done loading spectral cube.")
    self.save()

  def save(self):
      self.speak('attempting to save the cube of loaded, shifted, compiled spectra to {0}'.format(self.filename))
      tosave = {}
      for s in self.savable:
          tosave[s] = self.__dict__[s]
          self.speak('  including [{0}] in the saved cube structure'.format(s))
      np.save(self.filename, tosave)

  def load(self):
      self.speak('attempting to load previously saved cubes from...')
      self.speak('{0}'.format(self.filename))
      loaded = np.load(self.filename)[()]
      for s in self.savable:
          self.__dict__[s] = astropy.table.Table(loaded[s])
          self.speak('  loading [{0}] from the saved cube structure'.format(s))
      self.numberofstars = len(self.cubes)
      self.numberofwavelengths = len(self.spectral)
      self.numberoftimes = len(self.temporal)

  def movieCube(self, fps=30, bitrate=1800*20, stride=10):
      '''Create movie of the spectral cube.'''

      self.populate()
      metadata = dict(artist='Z.K.B.-T.')
      self.writer = matplotlib.animation.FFMpegWriter(fps=fps, metadata=metadata, bitrate=bitrate)
      filename = self.obs.extractionDirectory + 'cube_{0}stars_{1}spectra_{2}stride.mp4'.format(self.numberofstars, self.numberoftimes, stride)
      self.plotSpectra(0)
      with self.writer.saving(self.figure, filename, self.figure.get_dpi()):
          # loop over exposures
          for i in np.arange(self.numberoftimes)[::stride]:
              self.speak('plotting {0}/{1}'.format(i, self.numberoftimes))
              self.plotSpectra(i)
              self.writer.grab_frame()
      self.speak('saved movie to {0}'.format(filename))
      os.system('open {0}'.format(filename))

  def plotSpectra(self, which):
        '''For the ith component in the cube, plot all the spectra.'''
        self.populate()
        linekeys = ['airmass', 'rotatore', 'cosmics', 'sky', 'width', 'centroid', 'shift', 'lc']
        # set up a plot window to show how the interpolation is going
        self.figure = plt.figure('spectra', figsize=(24,12), dpi=100)
        # define an empty cube that we're going to populate with spectra
        try:
            self.ax_spectra
        except:
            self.ax_spectra, self.ps = {}, {}
            sharex=None
            gs = plt.matplotlib.gridspec.GridSpec(len(self.cubekeys),self.numberofstars,hspace=0.05,wspace=0.2, left=0.1,right=0.7)
            gstemporal = plt.matplotlib.gridspec.GridSpec(1, len(linekeys), left=0.75, right=0.95, wspace=0.05)
            sharex=None
            ok = self.temporal['ok']
            for s in range(self.numberofstars):
                fine = self.cubes['ok'][s,ok,:]
                for i in range(len(self.cubekeys)):
                    k = self.cubekeys[i]
                    if s ==0:
                        self.ax_spectra[k] = []
                    self.ax_spectra[k].append(plt.subplot(gs[i,s], sharex=sharex))
                    sharex = self.ax_spectra[k][0]
                    self.ax_spectra[k][0].set_ylabel(k)
                    self.ax_spectra[k][s].set_ylim(np.nanmin(self.cubes[k][s,ok,:][fine]),np.nanmax(self.cubes[k][s,ok,:][fine]) )
                    plt.setp(self.ax_spectra[k][s].get_xticklabels(), visible=False)

                plt.setp(self.ax_spectra[k][s].get_xticklabels(), visible=True)
            self.ax_spectra[k][0].set_xlim(np.min(self.spectral['wavelength']), np.max(self.spectral['wavelength']))

            kw = dict(color='black', linewidth=1)

            sharey = None
            self.bars = []
            for i in range(len(linekeys)):
                l = linekeys[i]
                ax = plt.subplot(gstemporal[i], sharey=sharey)
                sharey=ax
                ax.plot(self.temporal[l], np.arange(self.numberoftimes), alpha=0.25, **kw)
                ax.plot(self.temporal[l][ok], np.arange(self.numberoftimes)[ok], **kw)
                ax.set_xlabel(l,rotation=45)
                ax.set_xlim(np.nanmin(self.temporal[l][ok]),np.nanmax(self.temporal[l][ok]))
                self.bars.append(ax.axhline(which, color='gray', linewidth=4, alpha=0.25))
                plt.setp(ax.get_yticklabels(), visible=False)
                plt.setp(ax.get_xticklabels(), visible=False)
                if i == len(linekeys)/2:
                    self.timestamp = ax.set_title('{0}: {1} {2}'.format(self.obs.name, self.temporal['ut-date'][which], self.temporal['ut-time'][which]))
            ax.set_ylim(self.numberoftimes, -1)
            for s in range(self.numberofstars):
                for k in self.cubekeys:
                    if s ==0:
                        self.ps[k] = []
                    self.ps[k].append(self.ax_spectra[k][s].plot(self.spectral['wavelength'], self.cubes[k][s,which,:], **kw)[0])

        for b in self.bars:
            b.set_ydata(which)
        for s in range(self.numberofstars):
            for k in self.cubekeys:
                fine = self.cubes['ok'][s,which,:]
                self.ps[k][s].set_data(self.spectral['wavelength'][fine], self.cubes[k][s,which,:][fine])
                self.ps[k][s].set_alpha(self.temporal['ok'][which]*0.75 + 0.25)
            self.ax_spectra[self.cubekeys[0]][s].set_title('image {0}, star {1}'.format(self.temporal['n'][which], s))
        self.timestamp.set_text('{0}: {1} {2}'.format(self.obs.name, self.temporal['ut-date'][which], self.temporal['ut-time'][which]))

  def shiftCube(self,plot=False):
    '''Shift all the spectra for a star to line up in wavelength.'''
    self.speak('shifting all spectra to line up at the calcium triplet')


    if plot:
        plt.figure('shifting spectra')
        plt.ion()

    # select a narrow wavelength range near the Ca triplet
    left = np.argmin(np.abs(self.spectral['wavelength'] - 8350))
    right = np.argmin(np.abs(self.spectral['wavelength'] - 8800))
    width = 5
    triplet = [8498.,8542.,8662.]

    c = self.cubes['raw_counts']

    # create a template to shift relative to
    wave = self.spectral['wavelength'][left:right]
    start = np.ones_like(wave)*np.std(zachopy.oned.subtractContinuum(c[self.goodComps[0],c.shape[1]/2,left:right]))*3
    for t in triplet:
      start *= 1 - np.exp(-(wave - t)**2/2/width**2)
    start -= np.median(start)

    shifts = np.zeros((self.numberofstars, self.numberoftimes))
    for i in range(len(c[0,:,0])):
      for star in range(len(c[:,0,0])):
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
          for key in self.cubes.keys():
              interpolation = scipy.interpolate.interp1d(pixels,self.cubes[key][star,i,:], kind='linear', bounds_error=False, fill_value=0.0)
              self.cubes[key][star,i,:] = interpolation(pixels + offset)
              self.speak('  shifting [{0}] by {1}A'.format(key, offset))
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
            new = zachopy.oned.subtractContinuum(self.shifted['raw_counts'][star,i,left:right])
            ax[1].plot(wave, new, color='green', alpha=0.9, linewidth=2)
            ax[1].set_xlim(wave.min(), wave.max())
            ax[1].set_autoscaley_on
            ax[0].set_title('{0}/{1} stars; {2}/{3} spectra'.format(star+1, len(c[:,0,0]), i+1, len(c[0,:,0])))
            ax[0].axvline(len(x)/2.0)
            plt.draw()

          shifts[star,i] = offset
          self.speak( "shift = {4} for {0}/{1} stars; {2}/{3} spectra".format(star+1, len(c[:,0,0]), i+1, len(c[0,:,0]), offset))
    self.squares['shift'] = shifts
    self.temporal['shift'] = self.squares['shift'][self.obs.target[0],:]

  def convolveCube(self,width=0.5):
    '''Take a spectral cube, convolve it with a Gaussian (NOT USED!).'''
    c = self.flux
    w = self.spectral['wavelength']
    plt.ion()
    x = np.arange(-width*10, width*10)
    gauss = lambda x: np.exp(-0.5*x**2/width**2)
    new_cube = np.copy(c)
    for star in range(len(c[:,0,0])):
        plt.figure('convolution')

        for i in range(len(c[0,:,0])):
            new_cube[star,i,:] = np.convolve(c[star,i,:], gauss(x),'same')
            plt.plot(c[star,i,:])
            plt.plot(new_cube[star,i,:])

    self.flux = new_cube

  def createBins(self, binsize=250, remake=False):
      '''Take a spectral cube, and creates big wavelength bins.'''
      self.populate()
      self.speak("Binning the spectral cube to {0}A binsize.".format(binsize))
      #self.shiftCube(plot=True)
      #self.show()
      wavelength = self.spectral['wavelength']
      dnpixelsdw = self.spectral['dnpixelsdw']

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
      for k in self.cubes.keys():
          self.speak('binning {0} to {1}A bins'.format(k, binsize))
          shape = (self.numberofstars, self.numberoftimes, self.numberofwavelengths/binsize, binsize)

          # integrate spectrum over bin (dndw factor accounts for the fractional pixels represented by resampled cube spectrum
          # !!!!! DO I NEED TO INCLUDE DNPIXELSDW??
          # this was the old version:
          #                      self.binned_cubes[newkey][star, time, bin] = scipy.integrate.trapz(wholespectrum[mask]*dnpixelsdw[mask], wavelength[mask])
          if k == 'ok':
              dndw = 1
          else:
              dndw = self.spectral['dnpixelsdw'].reshape(1, 1, self.numberofwavelengths)
          self.binned_cubes[k] = (self.cubes[k]*dndw).reshape(shape).sum(-1)
          if k=='width' or k =='centroid' or k =='peak':
              self.binned_cubes[k] /= (dndw).reshape((1,1, self.numberofwavelengths/binsize, binsize)).sum(-1)
          if k=='ok':
              self.binned_cubes[k] = (self.binned_cubes[k]/binsize).astype(np.bool)
      #
      self.bin_centers = bin_centers
      self.binned_cubes = astropy.table.Table(self.binned_cubes)
      #np.save(binned_filename, (self.binned_cubes, self.bin_centers, self.binned_cubes['ok']))
      #print "    Saved binned cube to {0}".format(binned_filename)

  def correctBins(self):
    '''Use comparison stars to correct for atmospheric losses, create a self.binned_corrected.'''

    self.populate()
    self.speak('using comparison stars {0} to correct for atmospheric losses'.format(self.obs.goodComps))
    wavelength = self.bin_centers
    vmin = 0.98
    vmax = 1.02
    nStars, nTimes, nWaves = self.binned_cubes['raw_counts'].shape

    # create empty correction and uncertainty arrays
    correction, uncertainty = np.ones((nTimes, nWaves)), np.ones((nTimes, nWaves))

    def weightedsum(array):
        return ((array*self.binned_cubes['ok'])[self.goodComps,:,:].sum(0)/(self.binned_cubes['ok'])[self.goodComps,:,:].sum(0))

    correction = weightedsum(self.binned_cubes['raw_counts'])
    uncertainty = np.sqrt(weightedsum(self.binned_cubes['raw_counts'] + self.binned_cubes['sky']))

    self.binned_correction = np.ma.MaskedArray(correction, mask=((self.binned_cubes['ok']==False).sum(0).astype(np.bool)), fill_value=np.nan)
    self.binned_correction_uncertainty = np.ma.MaskedArray(uncertainty, mask=((self.binned_cubes['ok']==False).sum(0).astype(np.bool)), fill_value=np.nan)

    # normalize the correction spectrum to be close to one
    mediancompositespectrum = np.ma.median(self.binned_correction, 0)
    self.binned_correction /= mediancompositespectrum.reshape(1,nWaves)
    self.binned_correction_uncertainty /= mediancompositespectrum.reshape(1,nWaves)
    #self.display.one(self.binned_correction.filled(), clobber=True)
    #self.display.one(self.binned_correction_uncertainty.filled())
    self.binned_cubes['corrected'] = self.binned_cubes['raw_counts']/self.binned_correction
    photonnoise = np.sqrt(self.binned_cubes['raw_counts'] + self.binned_cubes['sky'])/self.binned_cubes['raw_counts']
    correctionnoise = self.binned_correction_uncertainty
    self.binned_cubes['uncertainty'] = np.sqrt(photonnoise**2 + correctionnoise**2)

  def makeMeanSpectrum(self, plot=False):
    #self.loadSpectra()
    self.populate()
    wavelength = self.spectral['wavelength']
    spectrum = np.median(self.cubes['raw_counts'][self.obs.target,:,:],1).flatten()
    assert(len(spectrum) == len(wavelength))
    if plot:
        fi, ax = plt.subplots(1)
        unit = 10.0
        ax.plot(wavelength/unit, spectrum*unit)
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Flux (photons/nm/exposure)')
    self.speak("saving median spectrum to")
    self.speak( self.obs.extractionDirectory + 'medianSpectrum.npy')
    np.save(self.obs.extractionDirectory + 'medianSpectrum.npy', (wavelength, spectrum))

  def makeLCs(self,binsize=250, remake=False):
    '''Wrapper to go from extracted spectra to binned, multiwavelength lightcurves.'''

    self.populate()
    # make (and save) and mean spectrum, for plotting comparisons at later steps
    self.makeMeanSpectrum()

    # pick the target star
    target = self.obs.target[0]
    comparisons = self.obs.goodComps

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



    lcDirectory = self.obs.extractionDirectory + 'chromatic' + ('%d' % binsize) + '/'
    zachopy.utils.mkdir(lcDirectory)
    lcDirectory = lcDirectory + 'originalLCs/'
    zachopy.utils.mkdir(lcDirectory)

    self.lcs = []

    # loop through wavelength bins
    for wave in range(nWaves):
      left, right = bin_starts[wave], bin_ends[wave]
      lcfilename =  lcDirectory + '/{0:05d}to{1:05d}.lightcurve'.format(np.int(left), np.int(right))

      # is there *any* good data at this wavelength?
      if self.binned_cubes['ok'][target,:,wave].any():

        # determine where the target star is ok and the correction is ok (on a time-by-time basis)
        ok = ((self.binned_cubes['ok'][target,:,wave] != False).flatten()*(self.binned_correction.mask[:,wave] == False).flatten())

        # report what's being made
        self.speak('making light curve for {0} to {1}A, with {2} good points'.format(bin_starts[wave], bin_ends[wave], np.sum(ok != False)))

        # create an empty LC object
        lc = astropy.table.Table()
        lc['bjd'] = self.temporal['bjd']
        lc['flux'] = self.binned_cubes['corrected'][target,:,wave].flatten()/np.median(self.binned_cubes['raw_counts'][target,:,wave].flatten())
        lc['uncertainty'] = self.binned_cubes['uncertainty'][target,:,wave].flatten()
        lc['ok'] = ok.astype(np.int)

        # pull out global values
        for key in ['airmass', 'rotatore']:
            lc['{0}'.format(key)] = self.temporal[key]

        # pull out the star-by-star (wavelength-independent) quantities
        for key in ['width', 'centroid', 'shift']:
            lc['{0}_target'.format(key)] = self.squares[key][target]
            for comparison in comparisons:
                lc['{0}_star{1:02.0f}'.format(key, comparison)] = self.squares[key][comparison]

        # pull out the star-by-star wavelength specific values
        for key in ['sky', 'peak']:
            lc['{0}_target'.format(key)] = self.binned_cubes[key][target,:,wave]
            for comparison in comparisons:
                lc['{0}_star{1:02.0f}'.format(key, comparison)] = self.binned_cubes[key][comparison,:,wave]

        # pull out the star-by-star wavelength specific values that should be measured relative to the more global values
        for key in ['width', 'centroid']:
            lc['d{0}_target'.format(key)] = self.binned_cubes[key][target,:,wave] - lc['{0}_target'.format(key)]
            for comparison in comparisons:
                lc['d{0}_star{1:02.0f}'.format(key, comparison)] = self.binned_cubes[key][comparison,:,wave] - lc['{0}_star{1:02.0f}'.format(key, comparison)]

        #lc.populate(bjd, flux, uncertainty, **lc)
        table = astropy.table.Table(lc)
        table['bjd'].format = '.10f'
        # MAYBE THIS IS A KLUDGE, MAYBE IT ISN'T?
        table = table[table['ok'].astype(np.bool)]

        table.write(lcfilename, format='ascii.fixed_width')
        self.speak('saved light curve to')
        self.speak('{0}'.format(lcfilename))

        '''for key in self.binned_cubes.keys():
        if key != 'corrected' and key != 'error':
          dict[key+'_target'] = self.binned_cubes[key][target,ok,wave].flatten()
          for comparison in self.obs.goodComps:
            dict[key+'_comparison{0:02.0f}'.format(comparison)] = self.binned_cubes[key][comparison,ok,wave].flatten()
        for k in keystoinclude:
          if k == 'ok':
              continue
          if k == 'airmass':
              newkey = k
          else:
              newkey = 'global_{0}'.format(k)
          dict[k] = self.temporal[k][ok]

          assert(np.isfinite(flux[ok]).all())
          assert(np.sum(ok) > 0)
          print bjd.flatten()[ok].size
          lc.populate(bjd[ok], flux[ok], error[ok], **dict)

          #lc.plot()
          print lc
          lc.save()
          self.lcs.append(lc)'''

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




  def imageTarget(self, title=None, vmin=None, vmax=None):
    '''Make an image of the input cube.'''
    if title is None:
      title = self.obs.name + ' | ' + self.obs.night
    print "   Trying to image bins for " + title
    bin_centers = self.bin_centers
    bin_ok = self.binned_cubes['ok']
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
    self.lc = np.zeros(bjd.size, types)
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
      self.cubekeystoplot = ['flux', 'raw_counts_target', 'sky_target', 'airmass', 'width_target', 'centroid_target', 'peak_target']
      gs = plt.matplotlib.gridspec.GridSpec(len(self.cubekeystoplot), 1, height_ratios=[5,1,1,1,1,1,1], wspace=0, hspace=0)
      self.ax = {}
      for i in range(len(self.cubekeystoplot)):
        key = self.cubekeystoplot[i]
        try:
          sharex = self.ax[self.cubekeystoplot[0]]
        except:
          sharex = None
        self.ax[key] = plt.subplot(gs[i], sharex=sharex)
    for key in self.cubekeystoplot:
      self.ax[key].plot(self.lc['bjd'], self.lc[key], markersize=1, marker='o', alpha=0.3, color='black', linewidth=0)
      self.ax[key].set_xlim(self.lc['bjd'].min(), self.lc['bjd'].max())
      self.ax[key].set_ylabel(key)

    plt.draw()
    bla = raw_input(self.filename + '?')
