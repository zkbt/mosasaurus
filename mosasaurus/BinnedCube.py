'''
This is still a sketch.
The goal is to to make this inherit from Cube,
but it's created via

c = Cube()
bc = BinnedCube(c, binsize=50)

but then has access to the same visualization, saving, and processing tools.
'''


'''
The stuff below is mostly code that came out of the original Cube (some not)


  ### FIX ME ### (looks like this was made before multiple widths became a thing)
  def roughLC(self, target=None, comps=None, wavelengths=None, **kwargs):
      '''construct a rough LC, over a given wavelength bin'''

      # if comps is just one-element, make it into a list
      target = self.target
      comps = self.comparisons

      self.speak('updating the rough LC with target {0} and comps {1} and wavelength range {2}'.format(target, comps, wavelengths))

      w = self.spectral['wavelength']
      if wavelengths is None:
          wavelengths = [np.min(w), np.max(w)]
      blueenough = (w <= np.max(wavelengths))
      redenough = (w >= np.min(wavelengths))
      waveok = blueenough*redenough
      targetflux = self.cubes['raw_counts'][target][:,waveok].sum(-1)
      comparisonflux = self.cubes['raw_counts'][comps].sum(0)[:,waveok].sum(-1)

      self.speak('comp/targ is typically {0}'.format(np.median(comparisonflux/targetflux)))
      self.temporal['lc'] = targetflux/comparisonflux
      self.temporal['lc'] /= np.median(self.temporal['lc'])
      self.speak('rms is {0}'.format(np.std(self.temporal['lc'])))

      plt.ion()
      x, y = self.temporal['bjd'], self.temporal['lc']
      ok = self.temporal['ok']
      plt.plot(x[ok], y[ok], **kwargs)

  # ????????
  def doubleplot(self,
                    binsize=500,
                    wavemin=4500, wavemax=9000,
                    divide=False, ylim=None,
                    median=False,
                    title=None):

        self.populate()
        c = self
        name = c.obs.target.name
        date = c.obs.night.name
        wavelengths = np.arange(wavemin, wavemax, binsize)



        lcs = []

        import craftroom.cmaps
        blue, red = 'indigo', 'darkorange'
        cmap = craftroom.cmaps.one2another(blue, red)

        plt.ioff()
        for i, w in enumerate(wavelengths):
            c.roughLC(wavelengths=[w-binsize/2, w+binsize/2])
            ok = np.array(c.temporal['ok'])
            lc = c.temporal['bjd']+0.0, c.temporal['lc'] +0.0, ok
            lcs.append(lc)

        offset = int(c.temporal['bjd'][0])
        plt.figure(name, figsize=(8,6), dpi=70)
        if median:
            gs = plt.matplotlib.gridspec.GridSpec(3,1,
                        height_ratios=[1,1,.5], hspace=0.02)
        else:
            gs = plt.matplotlib.gridspec.GridSpec(2,1,
                        height_ratios=[1,.5], hspace=0.02)
        axlc = plt.subplot(gs[-1])

        plt.xlabel('BJD - {0}'.format(offset))
        plt.ylabel('Relative Flux')


        n, m = len(lc[2]), len(wavelengths)
        image, imageok = np.zeros((n,m)),  np.zeros((n,m))

        for i, lc in enumerate(lcs):
            fraction = i/(len(wavelengths) - 1.0)
            t = lc[0] - offset
            flux = lc[1]
            ok = lc[2]
            image[:,i] = flux + 0.0
            imageok[:,i] = ok + 0
            #print ok
            plt.plot(t[ok], flux[ok], alpha=binsize/500.0, linewidth=1,color=cmap(fraction))

        if ylim is None:
            valid = np.nonzero((imageok*np.isfinite(image)).flatten())[0]
            ylim = np.percentile(image.flatten()[valid], [1,99])

        plt.ylim(*ylim)
        plt.xlim(t.min(), t.max())


        if divide:
            image /= np.median(image, 1)[:,np.newaxis]

        axim = plt.subplot(gs[0])
        kw = dict(interpolation='nearest', cmap='gray',
                                vmin=ylim[0], vmax=ylim[1], aspect='auto',
                                extent=[min(t), max(t),
                                        (min(wavelengths) - binsize/2.0)/10, (max(wavelengths) + binsize/2.0)/10],
                                origin='lower')
        axim.imshow(image.T, **kw)


        plt.setp(axim.get_xticklabels(), visible=False)
        if title is None:
            title = '{name} with {instrument}\n[from {0}nm ({blue}) to {1}nm ({red}) in {2}nm-wide bins]'.format( wavemin/10, wavemax/10, binsize/10, name=name,blue=blue, red=red, instrument=c.obs.instrument.name)
        plt.title(title)
        plt.ylabel('Wavelength (nm)')

        if median:
            axmed = plt.subplot(gs[1])
            divided = (image/np.median(image, 1)[:,np.newaxis])
            kw['vmin'], kw['vmax'] = np.percentile(divided, [5,95])
            axmed.imshow(divided.T, **kw)
            plt.setp(axmed.get_xticklabels(), visible=False)

        plt.draw()

        if divide:
            filename = '{0}_binto{1}_{2}_normalized.pdf'.format(name,binsize, date)
        else:
            filename = '{0}_binto{1}_{2}.pdf'.format(name,binsize, date)
        plt.savefig(filename)


  ### FIX ME (not used)
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
      #dnpixelsdw = self.spectral['dnpixelsdw']

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
          #if k == 'ok':
          #      dndw = 1
          #else:
          #  dndw = self.spectral['dnpixelsdw'].reshape(1, 1, self.numberofwavelengths)
          self.binned_cubes[k] = (self.cubes[k]*dndw).reshape(shape).sum(-1)
          if k=='width' or k =='centroid' or k =='peak':
              self.binned_cubes[k] /= (dndw).reshape((1,1, self.numberofwavelengths/binsize, binsize)).sum(-1)
          if k=='ok':
              self.binned_cubes[k] = (self.binned_cubes[k]/binsize).astype(np.bool)
      #
      self.bin_centers = bin_centers
      self.binned_cubes = astropy.table.Table(self.binned_cubes)
      #np.save(binned_filename, (self.binned_cubes, self.bin_centers, self.binned_cubes['ok']))
      #self.speak("    Saved binned cube to {0}".format(binned_filename))

  def correctBins(self, **kw):
    '''Use comparison stars to correct for atmospheric losses, create a self.binned_corrected.'''

    self.populate()
    self.speak('using comparison stars {0} to correct for atmospheric losses'.format(self.comparisons))
    wavelength = self.bin_centers
    vmin = 0.98
    vmax = 1.02
    nStars, nTimes, nWaves = self.binned_cubes['raw_counts'].shape

    # create empty correction and uncertainty arrays
    correction, uncertainty = np.ones((nTimes, nWaves)), np.ones((nTimes, nWaves))

    def weightedsum(array):
        return ((array*self.binned_cubes['ok'])[self.comparisons,:,:].sum(0)/(self.binned_cubes['ok'])[self.comparisons,:,:].sum(0))

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

  def imageBins(self, **kw):
      self.createBins(**kw)
      self.correctBins(**kw)

      figure = plt.figure(figsize=(10,10), dpi=70)
      gs = plt.matplotlib.gridspec.GridSpec(self.numberofstars, 2)

      kw = dict(cmap='gray')
      ax=None
      for i in range(self.numberofstars):
          ax = plt.subplot(gs[i,0], sharex=ax, sharey=ax)
          ax.imshow(self.binned_cubes['raw_counts'][i], **kw)
          ax = plt.subplot(gs[i,1], sharex=ax, sharey=ax)
          ax.imshow(self.binned_cubes['corrected'][i], **kw)

  def makeMeanSpectrum(self, plot=False):
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
    filename = os.path.join(self.directory, 'medianSpectrum.npy')
    self.speak(filename)
    np.save(filename, (wavelength, spectrum))

  def makeLCs(self,binsize=250, remake=False):
    '''Wrapper to go from extracted spectra to binned, multiwavelength lightcurves.'''

    self.populate()
    # make (and save) and mean spectrum, for plotting comparisons at later steps
    self.makeMeanSpectrum()

    # pick the target star
    target = self.target
    comparisons = self.obs.comparisons

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


    lcDirectory = os.path.join(self.directory,  "chromatic{binsize:05.0f}/".format(binsize=binsize))
    mkdir(lcDirectory)
    lcDirectory = os.path.join(lcDirectory, 'originalLCs/')
    mkdir(lcDirectory)

    self.lcs = []

    # loop through wavelength bins
    for wave in range(nWaves):
      left, right = bin_starts[wave], bin_ends[wave]
      lcfilename =  os.path.join(lcDirectory, '/{0:05d}to{1:05d}.lightcurve'.format(np.int(left), np.int(right)))

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
            try:
                lc['{0}_target'.format(key)] = self.squares[key][target]
                for comparison in comparisons:
                    lc['{0}_star{1:02.0f}'.format(key, comparison)] = self.squares[key][comparison]
            except KeyError:
                self.speak("{} couldn't be found!".format(key))
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

        #table = table[table['ok'].astype(np.bool)]
        # REMOVED TO MAKE SURE MASKING IS EASIER AT LATER STEP


        table.write(lcfilename, format='ascii.fixed_width', bookend=False)
        self.speak('saved light curve to')
        self.speak('{0}'.format(lcfilename))

        '''for key in self.binned_cubes.keys():
        if key != 'corrected' and key != 'error':
          dict[key+'_target'] = self.binned_cubes[key][target,ok,wave].flatten()
          for comparison in self.obs.comparisons:
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
          self.speak(bjd.flatten()[ok].size)
          lc.populate(bjd[ok], flux[ok], error[ok], **dict)

          #lc.plot()
          self.speak(lc)
          lc.save()
          self.lcs.append(lc)'''

  def loadLCs(self, binsize=250):
    lcDirectory = os.path.join(self.directory, 'lc_binby' + ('%d' % binsize) + '/')
    g = glob.glob(os.path.join(lcDirectory, 'lc_*.npy'))
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
      title = self.obs.target.name + ' | ' + self.obs.night
    self.speak("   Trying to image bins for " + title)
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
    plt.savefig(os.path.join(self.obs.reducer.extractionDirectory, + 'divided_{0}.pdf'.format(self.obs.night)))
    a = raw_input('imaging  target')

  def help(self):
      self.speak(cubehelp)


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
    self.filename = os.path.join(self.directory, 'lc_binby' + ('%d' % self.binsize) + '/lc_{0:05d}to{1:05d}.npy'.format(np.int(self.left), np.int(self.right)))

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
