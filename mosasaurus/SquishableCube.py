from .Cube import Cube
from .imports import *


# set a cmap for stars
starcm = craftroom.cmaps.one2another('royalblue', 'sienna')

class SquishableCube(Talker):
    def __init__(self, filename):

        self.filename = filename
        Talker.__init__(self)

        self.load()
        self.directory = os.path.split(self.filename)[0]
        self.shift = 'shifted.npy' in self.filename
        self.binned = False
        self.markBad()

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

        dir = os.path.join(self.directory, 'imagedcubes')
        mkdir(dir)
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

      okfilename = os.path.join(self.directory, 'okgrid.npy')

      try:
        self.cubes['ok'] = np.load(okfilename)[()]
        self.speak('loaded OK-ness from {}'.format(okfilename))
        assert(False)
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


    def squish(self, binning=50):
        '''
        Bin a cube down to a squished cube, using the ok as weights.
        '''
        self.speak('binning the cube!')

        originalwavelength = self.spectral['wavelength']
        dw = np.median(np.diff(originalwavelength))
        left = self.spectral['wavelength'][::binning]
        right = left + dw*binning
        centers = 0.5*(left + right)
        n = np.int(len(originalwavelength)/binning)

        self.binned_cube = {}
        self.binned_spectral = {}
        self.binned_cube['raw_counts'] = {}
        self.binned_cube['ok'] = {}

        for star in self.stars:
            self.binned_cube['raw_counts'][star] = np.zeros((self.numberoftimes, np.int(self.numberofwavelengths/binning) ))
            for i in range(self.numberoftimes):

                flux = self.cubes['raw_counts'][star][i,:]
                ok = self.cubes['ok'][star][i,:]

                binned_weighted = fluxconservingresample(
                        originalwavelength, flux*ok, centers)

                binned_weights = fluxconservingresample(
                        originalwavelength, ok, centers)

                final = binned_weighted/binned_weights
                assert(np.isfinite(final).any())
                self.binned_cube['raw_counts'][star][i,:] = final

                self.binned_cube['raw_counts'][star][i,:] = final

                #self.binned_cube['raw_counts'][star].append(final)
                #self.speak('{}'.format(i))
        #self.binned_cube['raw_counts'][star] = np.array(np.array(self.binned_cube['raw_counts'][star]))
        self.binned_spectral['wavelength'] = centers
        self.binned = True
