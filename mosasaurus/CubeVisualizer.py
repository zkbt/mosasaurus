from .imports import *


class CubeVisualizer(Talker):
    '''
    A tool for visualizing a squishable cube.
    '''
    def __init__(self, cube=None):
        '''
        Initialize a visualizer, and link it to a Squishable or Binned cube.
        '''

        Talker.__init__(self)
        self.cube = cube

    @property
    def loupe(self):
        try:
            return self._loupe
        except AttributeError:
            self._loupe = loupe()
            return self._loupe

    def explore(self, image=None, key='raw_counts', star=None, vmin=None, vmax=None, **kw):
        '''
        Visually explore the cube, displaying flux vs. time + wavelength.

            image = (ntimepoints x nwavelengths) image to visualize,
                    by imshowing the 2D array, and showing slices of it along
                    both the horizontal and vertical


        '''

        color = self.cube.starcolor(self.cube.target)

        if image is None:
            # make sure a star is defined
            if star == None:
                star = self.cube.target

            # do different things for different keys
            if key in self.cube.cubekeys:
                z = self.cube.cubes[key][star]
                vmin, vmax = 0, np.nanpercentile(z, 99)
            if key == 'corrected':
                z = self.corrected()
                vmin, vmax = 0.979, 1.01
            if key == 'wavelengthed':
                z = self.wavelengthed()
                vmin, vmax = 0.72, 1.1
        else:
            z = image
            if vmin is None:
                vmin = np.nanpercentile(z, 1)
            if vmax is None:
                vmin = np.nanpercentile(z, 99)

        aspect_ratio = 6.0/4.0
        tiny = 0.1
        width_ratios = [1.0 - tiny, tiny]
        height_ratios = [tiny*aspect_ratio, 1.0 - tiny*aspect_ratio]
        space = 0.05
        hspace = space*aspect_ratio
        wspace = space
        figsize= np.array([aspect_ratio, 1.0])*8


        wavelength = self.cube.spectral['wavelength']
        times = self.cube.temporal['bjd']
        times -= np.min(times)
        try:
            ok = self.ok
        except AttributeError:
            ok = np.ones_like(z).astype(np.bool)
        self.loupe.setup(z, ok=ok, yaxis=wavelength, xaxis=times,
                              aspect='auto',
                              hspace=hspace, wspace=wspace,
                              width_ratio=width_ratios,
                              height_ratios=height_ratios,
                              figsize=figsize,
                              labelfontsize=None,
                              datacolor=color,
                              crosshaircolor=color,
                              left=0.1, bottom=0.1,
                              **kw)
        self.loupe.ax['2d'].set_xlabel('Time from Start of Observation (days)')
        self.loupe.ax['2d'].set_ylabel('Wavelength ($\AA$)')
        self.loupe.ax['slicey'].set_ylim(np.max(wavelength), np.min(wavelength))

        # make the light curves into scattered points
        for thing in ['slicex', 'slicex_bad']:
            self.loupe.plotted[thing].set_marker('.')
            self.loupe.plotted[thing].set_linewidth(0)
        self.loupe.plotted['slicex_bad'].set_marker('x')
        #self.loupe.ax['slicex'].yaxis.tick_right()
        #self.loupe.ax['slicex'].yaxis.set_label_position('right')

        # make the lines thicker
        linewidth=3
        for line in ['slicey', 'crossy', 'crossyextend', 'crossx', 'crossxextend']:
            self.loupe.plotted[line].set_linewidth(linewidth)

        self.loupe.set_limits(vmin, vmax)
        plt.draw()
        #self.loupe.run()


    def corrected(self, key='raw_counts'):
        # normalize along the wavelength axis
        star = self.cube.target
        self.speak('calculating [{}] for [{}], corrected by the mega-calibrator'.format(key, star))
        target = self.cube.cubes[key][star]
        # KLUDGE!!!
        comparison = self.cube.cubes[key][self.cube.comparisons[0]]
        z = target/comparison
        oned = np.nanmedian(z, 0)
        z = z/oned[np.newaxis,:]
        return z

    def wavelengthed(self, key='raw_counts', star=None):
        if star == None:
            star = self.cube.target

        # normalize along the wavelength axis
        self.speak('calculating [{}] for [{}], normalized by its median spectrum'.format(key, star))
        z = self.cube.cubes[key][star]
        oned = np.nanmedian(z, 0)
        return z/oned[np.newaxis,:]

    def overlayQuality(self, color='tomato'):
        '''
        Plot the cosmic rays on top of the plot.

        ### KLUDGE -- test this; it's not clear it's working yet!
        '''
        cmap = craftroom.cmaps.one2another(color, color, alphabottom=1.0, alphatop=0.0)
        a = self.loupe.ax['2d']

        try:
            ok = self.ok
        except AttributeError:
            ok = self.cube.cubes['ok'][self.cube.target]

        self.overlay = a.imshow(ok.T,
                                    cmap=cmap,
                                    extent=self.loupe.extent,
                                    interpolation='nearest',
                                    zorder=5,
                                    origin='lower',
                                    aspect=self.loupe.ax['2d'].get_aspect()
        )

    """
    def zapCosmics(self, remake=False):
        # the
        cosmicfilename = os.path.join(self.cube.binneddirectory, 'cosmics.npy')
        try:
            notcosmics = np.load(cosmicfilename)[()]
            assert(remake == False)
        except (IOError, AssertionError):
            self.speak('trying to zap cosmics')
            z = self.corrected()
            filtered = scipy.signal.medfilt(z, (5, 15))
            cosmics = z - filtered

            wavelengthstd = 1.48*np.nanmedian(np.abs(cosmics - np.nanmedian(cosmics, 0)[np.newaxis,:]), 0)
            #mad(cosmics, 0)
            normalized = cosmics/wavelengthstd[np.newaxis, :]
            notcosmics = np.abs(normalized) < 5
            np.save(cosmicfilename, notcosmics)
        self.ok = notcosmics
        for star in self.cube.stars:
            self.cube.binned_cubes['ok'][star] *= self.ok
    """



    def makeSliceMovies(self, keys=['wavelengthed', 'corrected', 'raw_counts'],
                         remake=False,
                         stride=1):
        '''
        Make movies slicing through wavelength.
        '''

        wavelength = self.cube.spectral['wavelength']
        z = self.cube.cubes['raw_counts'][self.cube.target]
        spectrum = np.nanmedian(z, 0)

        for key in keys:

            if key == 'raw_counts':
                axislabel = 'Photons/$\AA$'
            else:
                axislabel = 'Relative\nFlux'

            self.explore(key=key)
            self.loupe.moveCrosshair(y=np.min(wavelength), x=None)
            self.loupe.plotted['slicey'].set_data(spectrum, wavelength)
            self.loupe.ax['slicey'].set_xlim(0, np.nanpercentile(spectrum, 99)*1.1)
            plt.setp(self.loupe.ax['slicey'].get_xticklabels(), visible=False)
            self.loupe.ax['slicex'].set_ylabel(axislabel)

            plotfilename = os.path.join(self.cube.directory, '{}.pdf'.format(key))
            plt.savefig(plotfilename, dpi=1000)

            filename = os.path.join(self.cube.directory, '{}.mp4'.format(key))
            self.speak('saving movie to {}'.format(filename))
            self.loupe.movieSlice(direction='y', filename=filename, remake=remake, stride=stride)
