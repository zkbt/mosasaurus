from imports import *
from Observation import Observation
import astropy.table, astropy.time
import matplotlib.cm

plt.ion()

# set a cmap for stars
starcm = zachopy.cmaps.one2another('magenta', 'limegreen')

class Cube(Talker):
  '''Cube object stores a -wavelength-flux datacube.'''
  def __init__(self, obs, remake=False, max=None, shift=True, **kwargs):
    '''Initialize a data cube and populate it with data.'''
    Talker.__init__(self, line=200, **kwargs)

    if type(obs) == str:
        self.obs = Observation(obs, nods9=True)
    else:
        self.obs = obs

    self.tempfilename = self.obs.extractionDirectory + 'tempSpectralCube.npy'
    self.cubekeys = ['sky',  'centroid', 'width', 'peak', 'raw_counts']

    self.savable = ['cubes', 'squares', 'temporal', 'spectral', 'stellar']
    self.sindex = 0
    self.tindex = 1
    self.windex = 2
    self.goodComps = self.obs.goodComps
    self.shift = shift

  def selectWidths(self):
    '''show the movies for the apertures, and decide on widths'''
    self.speak('please pick which extraction widths to use')
    for d in self.starDirectories:
        pass

  def populate(self, remake=False, max=None, visualize=True):
    try:
        self.cubes
    except:
        try:
            assert(remake==False)
            self.load()
        except:
            self.loadSpectra(remake=remake, max=max, visualize=visualize)
    self.markBad()

  @property
  def display(self):
    try:
        return self._display
    except:
        self._display = zachods9('cube')
        return self._display

  @property
  def starDirectories(self):
      return glob.glob(self.obs.extractionDirectory + 'aperture_*')


  @property
  def filename(self):
    if self.shift:
        return self.obs.extractionDirectory + 'shifted_spectralCube_{0}stars_{1}spectra.npy'.format(self.numberofstars, self.numberoftimes)
    else:
        return self.obs.extractionDirectory + 'unshifted_spectralCube_{0}stars_{1}spectra.npy'.format(self.numberofstars, self.numberoftimes)

  @property
  def stars(self):
      return self.stellar['aperture']

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
    # 1d, star
    self.stellar = {}

    # update
    self.speak("Loading the spectral cube.")
    # define the directories that contain extracted stellar spectra
    #self.starDirectories =
    self.numberofstars = len(self.starDirectories)
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
    self.headers = np.load(self.obs.workingDirectory + 'headers.npy')[()]

    self.stellar['aperture'] = [x.split('/')[-1] for x in self.starDirectories]
    # loop over the spectra
    for timepoint in range(self.numberoftimes):
        # loop over all the stars
        for istar, star in enumerate(self.stars):

          # ccd number for this image
          ccdn =  self.obs.nScience[timepoint]

          # find the available spectrum
          spectrumFile = self.starDirectories[istar] + '/supersampled{0:04.0f}.npy'.format(ccdn)
          extractedFile = self.starDirectories[istar] + '/extracted{0:04.0f}.npy'.format(ccdn)

          self.speak('trying to load {0}'.format(spectrumFile))

          # load the extracted spectrum (or truncate the cubes at this point)
          try:
              supersampled = np.load(spectrumFile)[()]
              self.speak('loaded {0}'.format(spectrumFile))
              extracted = np.load(extractedFile)[()]
              self.speak('loaded {0}'.format(extractedFile))
              #print supersampled.keys()
              #print extracted.keys()
          except IOError:
              # if we've run out of spectra to load, then truncate
              truncate = True
              self.speak('truncating cube at {0}'.format(spectrumFile))

              break
          # make sure we know how many aperture widths we're dealing with (for this spectrum)
          widths = np.sort([np.float(x.split('_')[-1].split('px')[0])
                                    for x in supersampled.keys()
                                        if 'raw_counts' in x])

          # loop over the measurement types and populate the cubes
          for key in self.cubekeys + ['ok']:

            for w in widths:
                # figure out the combined key
                widthkey = '{:04.1f}px'.format(w)
                # make sure a cube exists for this key
                try:
                    self.cubes[key]
                except KeyError:
                    self.cubes[key] = {}
                # make sure a cube exists for this star
                try:
                    self.cubes[key][star]
                except KeyError:
                    self.cubes[key][star] = {}
                # make sure a cube exists for this extraction width
                try:
                    self.cubes[key][star][widthkey]
                except KeyError:
                    if key == 'ok':
                        self.cubes[key][star][widthkey] = np.ones((self.numberoftimes, len(supersampled['wavelength']))).astype(np.bool)
                    else:
                        self.cubes[key][star][widthkey] = np.zeros((self.numberoftimes, len(supersampled['wavelength'])))
                        self.speak('updating cubes[{key}][{star}][{widthkey}][{timepoint},:]'.format(**locals()))

                # populate with the supersampled spectrum
                if key != 'ok':
                    self.cubes[key][star][widthkey][timepoint,:] = supersampled[key + '_' + widthkey]

                #print '!!!'
                if 'raw_counts' in key:
                    print sum(self.cubes[key][star][widthkey][timepoint,:])
                    assert(sum(self.cubes[key][star][widthkey][timepoint,:])!=5500)
                    assert(sum(self.cubes[key][star][widthkey][timepoint,:])>0.0)


          # pull out data from the (unsupersampled) spectra to populate a square with dimensions self.numberofstars x self.numberoftimes
          for w in widths:
              widthkey = '{:04.1f}px'.format(w)
              for key in ['sky', 'width', 'centroid', 'cosmicdiagnostic']:

                  try:
                      self.squares[key]
                  except KeyError:
                      self.squares[key] = {}
                  try:
                      self.squares[key][star]
                  except KeyError:
                      self.squares[key][star] = {}
                  try:
                      self.squares[key][star][widthkey]
                  except KeyError:
                      self.squares[key][star][widthkey] = np.zeros(self.numberoftimes)

                  self.squares[key][star][widthkey][timepoint] = np.median(extracted[w][key])
                  self.speak('updating squares[{key}][{star}][{widthkey}][{timepoint}]'.format(**locals()))

        # if we've run out of spectra, then break out of the loop (with truncated cubes)
        if truncate:
            break

        self.speak('{0}/{1} spectra loaded into cube'.format(timepoint, self.numberoftimes))
        # if the spectra for all stars were successfully loaded, then
        try:
            self.temporal['n']
        except KeyError:
            self.temporal['n'] = []
        self.temporal['n'].append(ccdn)

    # make sure everything is truncated properly
    if truncate:
        self.speak("couldn't find all requested spectra, so truncated cube at a length of {0}".format(timepoint))
        for key in self.cubes.keys():
            self.cubes[key] = self.cubes[key][star][widthkey][0:timepoint,:]
        for key in self.squares.keys():
            self.squares[key] = self.squares[key][star][widthkey][0:timepoint]


    #self.cubes = astropy.table.Table(self.cubes)
    #self.squares = astropy.table.Table(self.squares)

    self.temporal = astropy.table.join(self.headers, astropy.table.Table(self.temporal))


    #self.temporal['cosmicdiagnostic'] = self.squares['cosmicdiagnostic'].sum(0)
    #self.temporal['width'] = self.squares['width'][self.obs.target[0],:]
    #self.temporal['centroid'] = self.squares['centroid'][self.obs.target[0],:]
    #self.temporal['sky'] = self.squares['sky'][self.obs.target[0],:]


    self.temporal['ok'] = np.ones(self.numberoftimes).astype(np.bool)#self.temporal['cosmicdiagnostic'] < self.obs.cosmicAbandon
    #self.cubes['ok'] *= self.temporal['ok'].reshape(1,self.numberoftimes,1)

    #assert(len(self.cubes[0][0]) == len(self.squares[0][0]))
    #assert(len(self.cubes[0][0]) == len(self.temporal))



    # define some useful arrays
    self.spectral['wavelength'] = supersampled['wavelength']
    self.spectral['dnpixelsdw'] = supersampled['dnpixelsdw']
    self.numberofwavelengths = len(self.spectral['wavelength'])

    if self.shift:
        self.shiftCube(plot=visualize)
    self.speak("Done loading spectral cube.")

    #self.roughLC()

    self.save()

  def markBad(self):
      satlimit = 150000
      faintlimit = 1

      for star in self.stars:
          widths = np.sort([np.float(x.split('_')[-1].split('px')[0])
                                  for x in self.cubes['raw_counts'][star].keys()])

          for w in widths:
              widthkey = '{:04.1f}px'.format(w)

              # mark wavelengths where the width is zero (or nearby)
              buffer = 10 # go this many pixels beyond borders
              undefined = self.cubes['width'][star][widthkey].sum(0) > 0
              bad = np.convolve(undefined, np.ones(buffer).astype(np.float)/buffer, mode='same')
              self.cubes['ok'][star][widthkey] *= (bad >= 0.9)

              # mark saturated values as not cool
              #self.cubes['ok'][star][widthkey] *= (self.cubes['peak'][star][widthkey] < satlimit)

              # mark things exceed the cosmic threshold as not cool
              cosmics = self.squares['cosmicdiagnostic'][star][widthkey]
              self.cubes['ok'][star][widthkey] *= cosmics[:,np.newaxis] < self.obs.cosmicAbandon
              # (DOUBLE CHECK THIS ISN'T A KLUDGE!)

              # mark
              self.cubes['centroid'][star][widthkey] -= np.median(self.cubes['centroid'][star][widthkey][self.cubes['ok'][star][widthkey]])

  def roughLC(self, target=None, comps=None, wavelengths=None, **kwargs):
      if target is None:
          target = self.obs.target[0]
      if comps is None:
          comps = self.obs.goodComps
      try:
          len(comps)
      except TypeError:
          comps = [comps]
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

  def save(self):
      self.speak('attempting to save the cube of loaded, shifted, compiled spectra to {0}'.format(self.filename))
      tosave = {}
      for thing in self.savable:
          tosave[thing] = self.__dict__[thing]
          self.speak('  including [{0}] in the saved cube structure'.format(thing))
      np.save(self.filename, tosave)

  def load(self):
      self.speak('attempting to load previously saved cubes from...')
      self.speak('{0}'.format(self.filename))
      loaded = np.load(self.filename)[()]
      for thing in self.savable:
          self.__dict__[thing] = loaded[thing]
          self.speak('  loading [{0}] from the saved cube structure'.format(thing))
      self.numberofstars = len(self.stellar['aperture'])
      self.numberofwavelengths = len(self.spectral)
      self.numberoftimes = len(self.temporal)

  def doubleplot(self,
                    binsize=500,
                    wavemin=4500, wavemax=9000,
                    divide=False, ylim=None,
                    median=False,
                    title=None):

        self.populate()
        c = self
        name = c.obs.name
        date = c.obs.night
        wavelengths = np.arange(wavemin, wavemax, binsize)



        lcs = []

        import zachopy.cmaps
        blue, red = 'indigo', 'darkorange'
        cmap = zachopy.cmaps.one2another(blue, red)

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
            title = '{name} with {instrument}\n[from {0}nm ({blue}) to {1}nm ({red}) in {2}nm-wide bins]'.format( wavemin/10, wavemax/10, binsize/10, name=name,blue=blue, red=red, instrument=c.obs.instrument)
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

  def movieCube(self,   fps=30, # how many frames per second
                        bitrate=1800*20, # bitrate (this seems to work well),
                        stride=10):# each frame will skip over this many timepoints):
      '''Create movie of the spectral cube.'''

      # set up a movie file, which will be populated with plots
      metadata = dict(artist='Z.K.B.-T.')
      self.writer = matplotlib.animation.FFMpegWriter(fps=fps, metadata=metadata, bitrate=bitrate)

      # set the filename based on whether cube was shifted
      if self.shift:
          filename = self.obs.extractionDirectory + 'shifted_cube_{0}stars_{1}spectra_{2}stride.mp4'.format(self.numberofstars, self.numberoftimes, stride)
      else:
          filename = self.obs.extractionDirectory + 'cube_{0}stars_{1}spectra_{2}stride.mp4'.format(self.numberofstars, self.numberoftimes, stride)

      # plot the first spectrum, to set up the plots
      plt.ioff()
      self.plotSpectra(0, remake=True)

      # make the movie
      with self.writer.saving(self.figure, filename, self.figure.get_dpi()):
          # loop over exposures
          for i in np.arange(self.numberoftimes)[::stride]:
              self.speak('plotting {0}/{1}'.format(i, self.numberoftimes))
              self.plotSpectra(i)
              self.writer.grab_frame()

      # finish and display
      self.speak('saved movie to {0}'.format(filename))
      os.system('open {0}'.format(filename))


  def starcolor(self, s):
    '''return a color for a particular star'''
    number = {s:i for i,s in enumerate(self.stellar['aperture'])}[s]
    return starcm(number/(self.numberofstars-1.0))

  def plotSpectra(self, which, remake=False, wavelengthscale=10.0):
        '''For the ith component in the cube, plot all the spectra.'''

        # make sure the cube has been populated
        self.populate()

        # these (global) values will be plotted along the right side of the plot
        self.globallinekeys = ['airmass', 'rotatore']

        # these (star-by-star) values will be plotted along the right side
        self.starlinekeys = ['cosmicdiagnostic', 'sky', 'width', 'centroid']#, 'shift']#, 'lc']

        # these are the combination of linekeys
        self.linekeys = []
        self.linekeys.extend(self.starlinekeys)
        self.linekeys.extend(self.globallinekeys)

        # set up a plot figure
        self.figure = plt.figure('spectra', figsize=(8 + self.numberofstars*4,24), dpi=60)

        try:
            # determine whether or not the axes have been set up
            self.ax_spectra
            assert(remake == False)
        except (AssertionError, AttributeError):

            # create dictionaryies of axes and plots
            self.ax_spectra, self.ps = {}, {}

            # make gridspec structure with a row for each cubekey to plot, and a column for each star
            gs = plt.matplotlib.gridspec.GridSpec(len(self.cubekeys),self.numberofstars,hspace=0.08,wspace=0.2, left=0.04,right=0.75, bottom=0.05, top=0.95)

            # make a gridspec structure for vertical plots along the right side
            gstemporal = plt.matplotlib.gridspec.GridSpec(1, len(self.globallinekeys) + len(self.starlinekeys), left=0.78, right=0.98, wspace=0.05, bottom=0.05, top=0.95)
            sharex=None
            #ok = self.temporal['ok']

            # loop over the stars (columns)
            for istar, star in enumerate(self.stars):

                # figure out the widths
                widths = np.sort([np.float(x.split('_')[-1].split('px')[0])
                                        for x in self.cubes['raw_counts'][star].keys()])

                for i in range(len(self.cubekeys)):

                    # what's the key for the xaxis
                    k = self.cubekeys[i]


                    # add an axes for this key (row) and star (column)
                    try:
                        self.ax_spectra[k].append(plt.subplot(gs[i,istar], sharex=sharex))
                    except:
                        self.ax_spectra[k] = [plt.subplot(gs[i,istar], sharex=sharex)]
                    sharex = self.ax_spectra[k][0]
                    self.ax_spectra[k][0].set_ylabel(k)

                    # loop over widths, and plot them
                    for w in self.cubes['ok'][star].keys():

                        # a mask setting which data points are okay
                        fine = self.cubes['ok'][star][w][which,:]

                        # pull out the y-values for this spectrum plot
                        thisy = self.cubes[k][star][w][self.cubes['ok'][star][w]]
                        thisy = thisy[np.isfinite(thisy)]
                        # make some adjustments, set ylim, and remove xlabels
                        if k in ['sky', 'peak', 'raw_counts']:
                            ylim = (0, np.percentile(thisy,99.9)*1.2)
                        elif k in ['centroid']:
                            ylim =  np.percentile(thisy,[1,99])
                        elif k in ['width']:
                            ylim = (0,self.obs.widest/2.0)
                        self.ax_spectra[k][istar].set_ylim(*ylim )
                        plt.setp(self.ax_spectra[k][istar].get_xticklabels(), visible=False)

                # turn xlabels back on (just for bottom row)
                plt.setp(self.ax_spectra[k][istar].get_xticklabels(), visible=True)
                self.ax_spectra[k][istar].set_xlabel('Wavelength (nm)')

            # set the xlimits of the spectral plots
            self.ax_spectra[k][0].set_xlim(np.min(self.spectral['wavelength']/wavelengthscale), np.max(self.spectral['wavelength']/wavelengthscale))

            # now, plot the vertical plots along the right
            kw = dict(linewidth=1, alpha=1.0)
            sharey = None
            self.bars = []

            # do the global line keys
            for i, l in enumerate(self.linekeys):
                # create the appropriate axis
                ax = plt.subplot(gstemporal[i], sharey=sharey)
                sharey=ax
                if l in self.globallinekeys:
                    x = self.temporal[l]
                    ax.plot(x, np.arange(self.numberoftimes), color='black', **kw)
                    #ax.set_xlim(np.nanmin(self.temporal[l]),np.nanmax(self.temporal[l]))
                else:

                    for s in self.stars:
                        widths = np.sort([np.float(x.split('_')[-1].split('px')[0])
                                                for x in self.cubes['raw_counts'][s].keys()])
                        w = '{:04.1f}px'.format(np.max(widths))
                        x = self.squares[l][s][w]
                        if l == 'cosmicdiagnostic':
                            x = np.log(x)
                        ax.plot(x, np.arange(self.numberoftimes),
                                    color=self.starcolor(s),
                                    **kw)
                            #ax.set_xlim(np.nanmin(self.squares[l][s][w]),np.nanmax(self.squares[l][s][w]))

                # tidy up the plot
                ax.set_xlabel(l,rotation=45)
                self.bars.append(ax.axhline(which, color='gray', linewidth=4, alpha=0.5))
                plt.setp(ax.get_yticklabels(), visible=False)
                plt.setp(ax.get_xticklabels(), visible=False)
                ax.get_xaxis().get_major_formatter().set_useOffset(False)

                if i == len(self.linekeys)/2:
                    self.timestamp = ax.set_title('{0}: {1} {2}'.format(self.obs.name, self.temporal['ut-date'][which], self.temporal['ut-time'][which]))

            ax.set_ylim(self.numberoftimes, -1)
            for istar, star in enumerate(self.stars):
                self.ps[star] = {}
                for k in self.cubekeys:
                    try:
                        self.ps[star][k]
                    except KeyError:
                        self.ps[star][k] = {}
                    for w in self.cubes[k][star].keys():
                        self.ps[star][k][w] = self.ax_spectra[k][istar].plot(self.spectral['wavelength'], self.cubes[k][star][w][which,:], color=self.starcolor(star), **kw)[0]

        # set the position of the bars
        for b in self.bars:
            b.set_ydata(which)

        # set the spectra data for the plots
        for istar, s in enumerate(self.stars):
            for k in self.cubekeys:
                for w in self.cubes[k][s].keys():
                    fine = self.cubes['ok'][s][w][which,:]
                    x = self.spectral['wavelength'][fine]/wavelengthscale
                    y = self.cubes[k][s][w][which,:][fine]
                    self.ps[s][k][w].set_data(x, y)
                    #self.speak('median {} value for star {} is {}'.format(k, s, np.median(y)))
            self.ax_spectra[self.cubekeys[0]][istar].set_title('image {0},\nstar {1}, aperture {2}'.format(self.temporal['n'][which], istar, s.replace('aperture_', '')))
        self.timestamp.set_text('{0}: {1} {2}'.format(self.obs.name, self.temporal['ut-date'][which], self.temporal['ut-time'][which]))
        #plt.draw()
        #self.input('huh?')
  def oldformat(self, widths):
      pass

  def shiftCube(self,plot=False):
    '''Shift all the spectra for a star to line up in wavelength.'''
    self.speak('shifting all spectra to line up at the calcium triplet')


    if plot:
        plt.figure('shifting spectra')
        plt.ion()

    # select a narrow wavelength range near the Ca triplet
    l, r = self.obs.correlationRange

    left = np.argmin(np.abs(self.spectral['wavelength'] - l))
    right = np.argmin(np.abs(self.spectral['wavelength'] - r))
    width = self.obs.correlationSmooth
    triplet = self.obs.correlationAnchors

    c = self.cubes['raw_counts']

    # create a template to shift relative to (some fake, broadened absorption lines)
    wave = self.spectral['wavelength'][left:right]
    start = np.ones_like(wave).astype(np.float)
    for t in triplet:
      start *= 1 - np.exp(-(wave - t)**2/2/width**2)
    start -= np.median(start)


    # create an a structure to store the shifts in
    shifts = {}
    for star in self.stars:
        shifts[star] = {}
        for w in c[star].keys():
            shifts[star][w] = np.zeros(self.numberoftimes)

    # calculate all the shifts
    for i in range(self.numberoftimes):
      for star in self.stars:
          for w in c[star].keys():
              # pull out the spectrum for one star,
              # at one extraction width, at one time point,
              # and trim it to a narrow range around the correlation anchors
              spectrum = c[star][w][i,left:right]

              # subtract its continuum
              this = zachopy.oned.subtractContinuum(spectrum)

              # cross correlation with the anchors
              xc = np.correlate(this,start,'same')

              # do a quadratic fit to estimate the peak
              x = np.arange(len(xc))
              coeff = np.polynomial.polynomial.polyfit(x[xc.argmax()-5:xc.argmax()+5], xc[xc.argmax()-5:xc.argmax()+5], 2)
              fit = np.polynomial.polynomial.Polynomial(coeff)
              der = fit.deriv()
              peak = der.roots()

              # calculate the offset from the peak
              offset = peak - len(start)/2 + 1

              # now, use interpolation to shift all the arrays
              pixels = np.arange(self.numberofwavelengths)
              for key in self.cubes.keys():

                  # should this be 'linear' or 'nearest'?
                  interpolation = scipy.interpolate.interp1d(
                                    pixels,self.cubes[key][star][w][i,:],
                                    kind='linear',
                                    bounds_error=False,
                                    fill_value=0.0)

                  self.cubes[key][star][w][i,:] = interpolation(pixels + offset)
                  self.speak('shifting [{0}] by {1}A'.format(key, offset))

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
                new = zachopy.oned.subtractContinuum(self.cubes[key][star,i,left:right])
                ax[1].plot(wave, new, color='green', alpha=0.9, linewidth=2)
                ax[1].set_xlim(wave.min(), wave.max())
                ax[1].set_autoscaley_on
                ax[0].set_title('{0}/{1} stars; {2}/{3} spectra'.format(star+1, len(c[:,0,0]), i+1, len(c[0,:,0])))
                ax[0].axvline(len(x)/2.0)
                plt.draw()

                shifts[star,i] = offset
          self.speak( "shift = {4} for {0}/{1} stars; {2}/{3} spectra".format(star+1, len(c[:,0,0]), i+1, len(c[0,:,0]), offset))
    self.squares['shift'] = shifts
    #self.temporal['shift'] = self.squares['shift'][self.obs.target[0],:]


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

  def correctBins(self, **kw):
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



    lcDirectory = self.obs.extractionDirectory + "chromatic{binsize:05.0f}/".format(binsize=binsize)
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

  def help(self):
      print '''
Independental measured variables:

=======================
stored in cube.temporal
=======================
airmass(time)
this is the overall average airmass of the field

rotatore(time)
this the instrument rotator angle -- it probably tells us something about the
instrument's changing illumination and flexure.

=======================
stored in cube.squares
=======================

centroid(star,time)
the centroid of the star in the cross-dispersion direction. this variable is the
median of the centroids across all wavelengths, representing an overall shift of
the star's position through the spectrograph's optics

width(star,time)
the width of the star in the cross-dispersion direction. this variable is the
median of the widths across all wavelengths, representing an overall change,
most likely due to seeing.

shift(star,time)
by how much did we have to shift the star in the wavelength direction, to
make its spectral features line up with a reference spectrum. as the zeropoint
of the wavelength calibration is set by the position of the star in the focal
plane, this tracks the motion of the star in the dispersion direction (and
probably other stuff too)

=======================
stored in cube.cubes
=======================

raw_counts(star,time,wavelength)
the flux from each star, integrated over the cross-dispersion direction and with
flux from the diffuse sky subtracted. this is what we use to make light curves!

sky(star,time,wavelength)
the extracted sky brightness, which was subtracted in an average sense, from the
the 1D stellar spectrum. if this value was estimated incorrectly, it could
influence the chromatic light curves.

dcentroid(star,time,wavelength)
the centroid of the star in the cross-dispersion direction, measured at a
particular wavelength, relative to "centroid" (see above). this represents
changes in the position on the detection caused by things like internal flexure
or atmospheric refraction.

dwidth(star,time,wavelength)
the width of the star in the cross-dispersion direction, measured at a
particular wavelength, relative to "width". flux correlations with dwidth
could be diagnostic of problems with the sky subtraction, or the extraction

peak(star,time,wavelength)
the brightness of the brightest pixel in the cross-dispersion direction,
at each wavelength, for each star. it will correlate strongly with seeing,
but is a slightly different tracer. correlations with "peak" that cannot be
explained by seeing alone might point to problems in the detector non-linearity
'''
