from .Cube import Cube
from .imports import *
from .CubeVisualizer import CubeVisualizer

binnedsavable = ['binned_cubes', 'squares', 'temporal', 'binned_spectral', 'stellar', 'target', 'comparisons']

# set a cmap for stars
starcm = craftroom.cmaps.one2another('royalblue', 'sienna')

class SquishableCube(Talker):
    def __init__(self, filename):

        self.filename = filename
        Talker.__init__(self)

        self.load()
        self.basedirectory = os.path.split(self.filename)[0]
        self.shift = 'shifted.npy' in self.filename
        self.binned = False
        self.binsize = None
        self.markBad()
        self.visualizer = CubeVisualizer(self)

    '''
    def bin(self, binning):

        # make sure things line up
        assert(self.numberofwavelengths % binning == 0)

        binned_spectral = {}
        for k in self.spectral.keys():
            binned_spectral[k] =
        newwavelength =
        self.unbinned
    '''

    @property
    def binlabel(self):
        return 'binby{}'.format(self.binsize)
    @property
    def binneddirectory(self):
        #if self.binned:
        d = os.path.join(self.basedirectory, self.binlabel)
        mkdir(d)
        return d
        #else:
        #    return self.basedirectory




    # keep - load the entire cube from a single file
    def load(self):
        self.speak('attempting to load previously saved cubes from...')
        self.speak('{0}'.format(self.filename))
        loaded = np.load(self.filename)[()]
        self.savable = loaded.keys()
        for thing in self.savable:
            self.__dict__[thing] = loaded[thing]
            self.speak('  loading [{0}] from the saved cube structure'.format(thing))
        self.numberofstars = len(self.stellar['aperture'])
        self.numberofwavelengths = len(self.spectral['wavelength'])
        self.numberoftimes = len(self.temporal)
        self.cubekeys = self.cubes.keys()
        self.stars = self.stellar['aperture']

    # keep -- gives each star a unique color
    def starcolor(self, s):
        '''return a color for a particular star'''
        number = {s:i for i,s in enumerate(self.stellar['aperture'])}[s]
        return starcm(number/(self.numberofstars-1.0))


    @property
    def medianSpectrum(self):
        return np.median(self.cubes['raw_count'][self.target])

    @property
    def cubeMegaComparison(self):
        '''
        Return a fake cube entry for an combined mega-comparison star.
        '''
        key = 'raw_counts'
        if len(self.comparisons) == 1:
            d = {}
            star = self.comparisons[0]
            d[key] = self.cubes[key][star]
            return d
        else:
            raise NameError("Darn it -- the mega-comparison hasn't been implemented yet!")

    def corrected(self, key='raw_counts'):
        # normalize along the wavelength axis
        star = self.target
        self.speak('calculating [{}] for [{}], corrected by the mega-calibrator'.format(key, star))
        target = self.cubes[key][star]
        comparison = self.cubeMegaComparison[key]
        z = target/comparison
        oned = np.median(z, 0)
        z = z/oned[np.newaxis,:]
        return z

    def wavelengthed(self, key='raw_counts', star=None):
        if star == None:
            star = self.target

        # normalize along the wavelength axis
        self.speak('calculating [{}] for [{}], normalized by its median spectrum'.format(key, star))
        z = self.cubes[key][star]
        oned = np.median(z, 0)
        return z/oned[np.newaxis,:]

    def timed(self, key='raw_counts', star=None):
        if star == None:
            star = self.target

        # normalize along the wavelength axis
        self.speak('calculating [{}] for [{}], normalized by its median time series'.format(key, star))
        z = self.cubes[key][star]
        oned = np.median(z, 1)
        return z/oned[:, np.newaxis]

    def wavelengthedtimed(self, key='raw_counts', star=None):
        if star == None:
            star = self.target

        # normalize along the wavelength axis
        self.speak('calculating [{}] for [{}], normalized by its median spectrum, and then by its median time series'.format(key, star))
        z = self.cubes[key][star]
        medianspectrum = np.median(z, 0)
        wnormalized = z/medianspectrum[np.newaxis,:]
        mediantimeseries = np.median(wnormalized, 1)

        return wnormalized/mediantimeseries[:, np.newaxis]


    # keep - makes tidy image plots of spectral quantities
    def imageCube(self, keys=None, stars=None, remake=False):
        '''
        Show an imshow of every cube key.
        '''

        if keys == None:
            keys = self.cubekeys

        if stars == None:
            stars = self.stars

        # we'll plot various types of images
        options = {
                'basic':'Raw Extracted Quantities',
                'wavelength':'Wavelength-normalized Quantities',
                'comparison':'Comparison-divided Quantities'
              }

        d = os.path.join(self.binneddirectory, 'imagedcubes')
        mkdir(d)
        for option, description in options.items():
            filename = os.path.join(dir, '{}+{}+{}+{}.pdf'.format(option, {True:'shifted', False:'raw'}[self.shift], '-'.join(keys), '+'.join(stars)))
            if os.path.exists(filename) and (remake == False):
                self.speak('{} already exists'.format(filename))
                #assert(False)
                continue
            nrows = len(keys)
            ncols = len(stars)

            fi, ax = plt.subplots(nrows, ncols,
                                    sharex=True, sharey=True, figsize=(12,8),
                                    gridspec_kw=dict(hspace=0.1, wspace=0.02))
            plt.suptitle('{}, [{}shifted]'.format(description, {True:'', False:'un'}[self.shift]))

            w = self.spectral['wavelength']
            t = np.arange(self.numberoftimes)

            # set up the imshow parameters
            imkw = dict(
                extent = [np.min(w), np.max(w), np.min(t), np.max(t)],
                cmap = 'gray',
                interpolation='nearest',
                aspect = 'auto',
                origin = 'lower'
            )

            # is this just a single panel?
            single = len(keys)*len(stars) == 1
            for i, key in enumerate(keys):
                for j, star in enumerate(stars):
                    if single:
                        a = ax
                    else:
                        a = ax[i,j]
                    a.cla()

                    if option == 'basic':
                        # don't modify the measurements at all
                        self.speak('displaying the raw measurements for {}'.format(key))
                        z = self.cubes[key][star]
                        vmin, vmax = np.min(z), np.max(z)

                    if option == 'wavelength':
                        # normalize along the wavelength axis
                        self.speak('displaying the measurements normalized by their median spectrum for {}'.format(key))
                        z = self.cubes[key][star]
                        oned = np.median(z, 0)
                        z = z/oned[np.newaxis,:]
                        vmin, vmax = 0.8, 1.2

                    if option == 'comparison':
                        # divide by the comparison star[s]

                        self.speak('displaying the measurements divided by the comparison[s] for {}'.format(key))
                        target = self.cubes[key][star]
                        comparison = self.cubeMegaComparison
                        z = target/comparison
                        oned = np.median(z, 0)
                        z = z/oned[np.newaxis,:]
                        vmin, vmax = 0.98, 1.02

                    #vmin, vmax = np.percentile(z, [1,99])
                    self.speak('the limits for {} on {} are [{} to {}]'.format(key, star, vmin, vmax))
                    #self.input('are these OK?')
                    a.imshow(z, vmin=vmin, vmax=vmax, **imkw)

                    # fuss with the axis labels
                    if j == 0:
                        a.set_ylabel('{}\n(timepoints)'.format(key))
                    else:
                        plt.setp(a.get_yticklabels(), visible=False)
                    if i == 0:
                        a.set_title('{}'.format(star))

                    if i == (len(self.cubekeys)-1):
                        a.set_xlabel('Wavelength (angstroms)')
                    else:
                        plt.setp(a.get_xticklabels(), visible=False)

            plt.draw()
            plt.savefig(filename, dpi=600)
            self.speak('saved image of this cube to {}'.format(filename))


    def markBad(self):
      '''mark bad time-wavelength-star data points as bad'''

      okfilename = os.path.join(self.basedirectory, 'okgrid.npy')

      try:
        self.cubes['ok'] = np.load(okfilename)[()]
        self.speak('loaded OK-ness from {}'.format(okfilename))
      except (IOError, AssertionError):

        self.cubes['ok'] = {}
        for star in self.stars:
            self.speak('marking bad points for {}'.format(star))
            self.speak('out of {} pixels...'.format(self.numberofwavelengths*self.numberoftimes))
            # what stars have no width defined for them?
            haswidth = self.cubes['width'][star].sum(0) > 0
            buffer = 10 # go this many wavelength pixels beyond borders
            smoother = np.convolve(haswidth, np.ones(buffer).astype(np.float)/buffer, mode='same')
            haswidth = (smoother >= 0.9)
            self.speak('   {} have no width'.format(np.sum(haswidth == False)))

            # what pixels are saturated (or the row is mostly saturated?)
            saturationlimit = 100000
            unsaturated = self.cubes['peak'][star] < saturationlimit
            self.speak('   {} are saturated'.format(np.sum(unsaturated == False)))


            # what times are weird outliers in this star's measurements?
            oned = np.ones_like(unsaturated).astype(np.bool)
            for key in ['centroid', 'width'] + self.shift*['shift']:
                x = self.squares[key][star]
                x -= np.median(x)
                oned *= (np.abs(x) < 10*1.48*mad(x))[:,np.newaxis]
            self.speak('   {} are weird outliers'.format(np.sum(oned == False)))

            # is finite
            finite = np.isfinite(self.cubes['raw_counts'][star])
            self.speak('   {} are infinite'.format(np.sum(finite == False)))


            # KLUDGE! this needs to handle each star individually!
            # isn't a cosmic ray (KLUDGE)
            '''
            if star == self.target:
                z = self.corrected()
                filtered = scipy.signal.medfilt(z, (3, 49))
                cosmics = z - filtered
                wavelengthstd = np.nanstd(cosmics, 0)
                normalized = cosmics/wavelengthstd[np.newaxis, :]
                notcosmics = np.abs(normalized) < 4.5
            self.speak('   {} are cosmic rays'.format(np.sum(notcosmics == False)))
            '''

            ok = haswidth*unsaturated*oned*finite#*notcosmics

            # store this value
            self.cubes['ok'][star] = ok



            pixelfraction = np.sum(self.cubes['ok'][star] > 0).astype(np.float)/self.numberofwavelengths/self.numberoftimes
            exposurefraction = np.sum(np.max(self.cubes['ok'][star], 1) > 0)/self.numberoftimes
            wavelengthfraction = np.sum(np.max(self.cubes['ok'][star], 0) > 0)/self.numberofwavelengths

            self.speak('For {}...,\n   {:.3%} of pixels are OK\n   {:.3%} of exposures are at least partially OK\n   {:.3%} of wavelengths are at least partially OK'.format(star, pixelfraction, exposurefraction, wavelengthfraction))
            self.speak('{} end ok'.format(np.sum(self.cubes['ok'][star])))

        np.save(okfilename, self.cubes['ok'])
        self.speak('loaded OK-ness from {}'.format(okfilename))


    def squish(self, binning=50, remake=False):
        '''
        Bin a cube down to a squished cube, using the ok as weights.
        '''
        self.speak('binning the cube!')
        self.binsize = binning
        self.binnedfilename = os.path.join(self.binneddirectory, self.binlabel + 'cube.npy')

        #try:
        #    assert(remake == False)
        #    self.binned_cubes, self.binned_spectral = np.load(self.binnedfilename)[()]
        #    self.speak('loaded binned cube from {}'.format(self.binnedfilename))
        #except (IOError, AssertionError):
        if True:

            originalwavelength = self.spectral['wavelength']
            dw = np.median(np.diff(originalwavelength))
            left = self.spectral['wavelength'][::binning]
            right = left + dw*binning
            centers = 0.5*(left + right)
            n = np.int(len(originalwavelength)/binning)

            self.binned_cubes = {}
            self.binned_spectral = {}
            self.binned_cubes['ok'] = {}
            for star in self.stars:
                self.binned_cubes['ok'][star] = np.ones((self.numberoftimes, np.int(self.numberofwavelengths/binning) ))

            for k in ['raw_counts', 'sky', 'centroid', 'width', 'peak']:
                self.binned_cubes[k] = {}
                self.speak(' (currently on {})'.format(k))
                for star in self.stars:
                    self.binned_cubes[k][star] = np.zeros((self.numberoftimes, np.int(self.numberofwavelengths/binning) ))
                    for i in tqdm(range(self.numberoftimes)):

                        # pull out the good elements
                        y = self.cubes[k][star][i,:]
                        ok = self.cubes['ok'][star][i,:].astype(np.float)

                        # make a weighted resample
                        binned_weighted = fluxconservingresample(
                                originalwavelength, y*ok, centers)

                        # what were the weights?
                        binned_weights = fluxconservingresample(
                                originalwavelength, ok, centers)

                        # what would the weights be if everything were perfect?
                        binned_ones= fluxconservingresample(
                                originalwavelength, np.ones_like(ok), centers)

                        # calculate the final binned
                        final = binned_weighted/binned_weights

                        # make sure that something is finite
                        assert(np.isfinite(final).any())

                        # store the final one
                        self.binned_cubes[k][star][i,:] = final
                        self.binned_cubes['ok'][star][i,:] = binned_weights/binned_ones#*= binned_weights > 0
                        #self.speak(repr(self.binned_cubes['ok'][star][i,:]))
                    assert(np.sum(self.binned_cubes['ok'][star]) > 0)
                    #self.speak(np.sum(self.binned_cubes['ok'][star]))
                    #self.binned_cubes['raw_counts'][star].append(final)
                    #self.speak('{}'.format(i))
            #self.binned_cubes['raw_counts'][star] = np.array(np.array(self.binned_cubes['raw_counts'][star]))
            self.binned_spectral['wavelength'] = centers
            self.binned_spectral['fractionofapixel'] = fluxconservingresample(
                    originalwavelength, self.spectral['fractionofapixel'], centers)


            for fraction in np.linspace(0, 1, 11):
                self.speak('{} points are <={} okay'.format(np.sum(self.binned_cubes['ok'][self.target] <= fraction), fraction))


            np.save(self.binnedfilename, (self.binned_cubes, self.binned_spectral))
            self.speak('saved binned cube to {}'.format(self.binnedfilename))
        self.binned = True


    def export(self):
        '''Export the data in this cube to an .npy.'''

        self.speak('attempting to save the binned cube of spectra to...')
        self.speak(' {0}'.format(self.binnedfilename))


        # create a dictionary, that's going to be saved
        tosave = {}
        tosave['cubes'] = self.binned_cubes
        tosave['spectral'] = self.binned_spectral
        othersavable = [ 'squares', 'temporal',  'stellar', 'target', 'comparisons']

        for thing in othersavable:
            tosave[thing] = self.__dict__[thing]

        for thing in tosave.keys():
            self.speak('  including [{0}] in the binned saved cube structure'.format(thing))

        # save that to a .npy file
        np.save(self.binnedfilename, tosave)


    def zapCosmics(self, wavelengthbin = 10, threshold = 5, plot=True, remake=False):
        '''
        Identify outliers in the divided cube, and call them cosmic rays.
        '''

        cosmicfilename = os.path.join(self.basedirectory, 'cosmicrays.npy')
        try:
            nocosmic = np.load(cosmicfilename)
            assert(remake == False)
        except (IOError, AssertionError):

            z = self.corrected()
            ntimes, nwavelengths = z.shape
            binned = np.mean(self.corrected().reshape(ntimes, np.ceil(nwavelengths/wavelengthbin).astype(np.int), wavelengthbin), 2)
            mediantimeseries = np.median(binned, 1)
            withouttimeseries = binned/mediantimeseries[:, np.newaxis] - 1
            wavelengthstd = np.nanstd(withouttimeseries, 0)
            normalized = withouttimeseries/wavelengthstd[np.newaxis, :]
            bad = np.abs(normalized) > threshold
            badoriginalshape = (bad[:,:,np.newaxis]*np.ones(wavelengthbin)).reshape((ntimes, nwavelengths))
            # kludge (individual stars should be different!)
            nocosmic = badoriginalshape == False
            if plot:
                filename = os.path.join(self.basedirectory, 'cosmicrays.pdf')
                self.visualizer.ok = badoriginalshape == False
                self.visualizer.explore(key='corrected', vmin=0.97, vmax=1.01, interpolation='nearest')
                self.visualizer.overlayQuality()
                plt.savefig(filename)
                self.speak('saved plot of cosmic identifications to {}'.format(filename))

            np.save(cosmicfilename, nocosmic)


        for star in self.stars:
            self.cubes['ok'][star] *= nocosmic


        '''
        self.squished = SquishedCube(   filename=self.binnedfilename,
                                        remake=True,
                                        cubes=self.binned_cubes,
                                        squares=self.squares,
                                        temporal=self.temporal,
                                        spectral=self.binned_spectral,
                                        stellar=self.stellar,
                                        target=self.target,
                                        comparisons=self.comparisons,
                                        directory=self.binneddirectory,
                                        binsize=self.binsize
                                        )
        return self.squished
        '''

class BinnedCube(SquishableCube):
    def __init__(self, filename=None,  **kwargs):
        '''
        This object stores a pre-written cube file.
        '''

        Talker.__init__(self)
        self.filename = filename
        self.load()
        self.binned = True
        self.binsize = np.int(filename.split('binby')[-1].split('cube')[0])
        self.directory = os.path.dirname(filename)
        self.visualizer = CubeVisualizer(self)


    def load(self):
        '''Load the data for this chromatic light curve from a .npy.'''

        self.speak('attempting to load binned cube from...')
        self.speak(' {0}'.format(self.filename))

        # load the dictionary
        loaded = np.load(self.filename)[()]
        self.savable = loaded.keys()
        for thing in self.savable:
            self.__dict__[thing] = loaded[thing]
            self.speak('  loading [{0}] from the binned cube file'.format(thing))

        self.numberofstars = len(self.stellar['aperture'])
        self.numberofwavelengths = len(self.spectral['wavelength'])
        self.numberoftimes = len(self.temporal)
        self.cubekeys = self.cubes.keys()
        self.stars = self.stellar['aperture']

    def makeLCs(self):
        '''Make a set of light curves from this binned cube.'''
        binsize = self.binsize
        lcdirectory = os.path.join(os.path.dirname(self.filename), 'lightcurves')
        mkdir(lcdirectory)
        # all of these are photons/angstron, so need to multiply by the binwidth
        # KLUDGE! this should be calculated for all stars, in the binning process!
        # calculate the noise for each light curve point!
        noisecountstarget = binsize*(self.cubes['raw_counts'][self.target] + self.cubes['sky'][self.target])
        signalcountstarget = binsize*(self.cubes['raw_counts'][self.target])
        fractionaltargetuncertainty = np.sqrt(noisecountstarget)/signalcountstarget

        noisecountscomparisons = np.zeros_like(noisecountstarget)
        signalcountscomparisons = np.zeros_like(noisecountstarget)
        for c in self.comparisons:
            noisecountscomparisons += binsize*(self.cubes['raw_counts'][c] + self.cubes['sky'][c])
            signalcountscomparisons += binsize*(self.cubes['raw_counts'][c])
        fractionalcomparisonsuncertainty = np.sqrt(noisecountscomparisons)/signalcountscomparisons

        uncertainty = np.sqrt(fractionaltargetuncertainty**2 + fractionalcomparisonsuncertainty**2)
        # this should be ntimes x nwavelengths shape


        wavelengths = self.spectral['wavelength']

        corrected = self.corrected()
        for i, w in enumerate(wavelengths):
            left = w - binsize/2
            right = w + binsize/2

            # create an empty LC object
            lc = astropy.table.Table()
            lc['bjd'] = self.temporal['bjd']

            lc['flux'] = corrected[:,i]
            lc['uncertainty'] = uncertainty[:,i]
            print("The median uncertainty is {}".format(np.median(lc['uncertainty'])))

            lc['ok'] = self.cubes['ok'][self.target][:,i]

            # pull out global time-dependent values
            for key in ['airmass', 'rotatore']:
                lc['{0}'.format(key)] = self.temporal[key]

            def starname(comparison):
                return comparison.replace('_', '-').replace('aperture', 'comp')

            # pull out the star-by-star (wavelength-independent) quantities
            for key in ['width', 'centroid', 'shift']:
                try:
                    lc['{0}_target'.format(key)] = self.squares[key][self.target]
                    for comparison in self.comparisons:
                        lc['{0}_{1}'.format(key, starname(comparison))] = self.squares[key][comparison]
                except KeyError:
                    self.speak("{} couldn't be found!".format(key))

            # pull out the star-by-star wavelength specific values
            for key in ['sky', 'peak']:
                lc['{0}_target'.format(key)] = self.cubes[key][self.target][:,i]
                for comparison in self.comparisons:
                    lc['{0}_{1}'.format(key, starname(comparison))] = self.cubes[key][comparison][:,i]

            # pull out the star-by-star wavelength specific values that should be measured relative to the more global values
            for key in ['width', 'centroid']:
                lc['d{0}_target'.format(key)] = self.cubes[key][self.target][:,i] - lc['{0}_target'.format(key)]
                for comparison in self.comparisons:
                    lc['d{0}_{1}'.format(key, starname(comparison))] = self.cubes[key][comparison][:,i] - lc['{0}_{1}'.format(key, starname(comparison))]

            #lc.populate(bjd, flux, uncertainty, **lc)
            table = astropy.table.Table(lc)
            table['bjd'].format = '.10f'

            #table = table[table['ok'].astype(np.bool)]
            # REMOVED TO MAKE SURE MASKING IS EASIER AT LATER STEP

            lcfilename =  os.path.join(lcdirectory, '{0:05d}to{1:05d}.lightcurve'.format(np.int(left), np.int(right)))

            table.write(lcfilename, format='ascii.fixed_width', bookend=False, overwrite=True)
            self.speak('saved light curve to')
            self.speak('{0}'.format(lcfilename))





    """
    def zapCosmics(self):
        # the
        z = self.corrected()
        filtered = scipy.signal.medfilt(z, (5, 15))
        cosmics = z - filtered

        wavelengthstd = 1.48*np.nanmedian(np.abs(cosmics - np.nanmedian(cosmics, 0)[np.newaxis,:]), 0)
        #mad(cosmics, 0)
        normalized = cosmics/wavelengthstd[np.newaxis, :]
        notcosmics = np.abs(normalized) < 5
        self.ok = notcosmics
        for star in self.squishable.stars:
            self.cubes['ok'][star] *= self.ok

    def corrected(self, key='raw_counts'):
        # normalize along the wavelength axis
        star = self.target
        self.speak('calculating [{}] for [{}], corrected by the mega-calibrator'.format(key, star))
        target = self.cubes[key][star]
        comparison = self.cubeMegaComparison[key]
        z = target/comparison
        oned = np.median(z, 0)
        z = z/oned[np.newaxis,:]
        return z

    @property
    def cubeMegaComparison(self):
        '''
        Return a fake cube entry for an combined mega-comparison star.
        '''
        key = 'raw_counts'
        if len(self.comparisons) == 1:
            d = {}
            star = self.comparisons[0]
            d[key] = self.cubes[key][star]
            return d
        else:
            raise NameError("Darn it -- the mega-comparison hasn't been implemented yet!")
    """
