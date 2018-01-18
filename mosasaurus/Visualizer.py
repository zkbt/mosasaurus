from .imports import *
from .SquishableCube import SquishableCube


class Visualizer(Talker):
    def __init__(self, squishablefilename):
        Talker.__init__(self)
        # load the cube of data
        self.squishable = SquishableCube(squishablefilename)

        #
        self.loupe = loupe()

    def explore(self, image=None, key='raw_counts', star=None, **kw):
        '''
        Set the image (key and star) to show.
        '''

        color = self.squishable.starcolor(self.squishable.target)

        if image is None:
            # make sure a star is defined
            if star == None:
                star = self.squishable.target

            # do different things for different keys
            if key in self.squishable.cubekeys:
                z = self.cubes[key][star]
                vmin, vmax = 0, np.nanpercentile(z, 99)
            if key == 'corrected':
                z = self.corrected()
                vmin, vmax = 0.979, 1.01
            if key == 'wavelengthed':
                z = self.wavelengthed()
                vmin, vmax = 0.72, 1.1
        else:
            z = image

        aspect_ratio = 6.0/4.0
        tiny = 0.1
        width_ratios = [1.0 - tiny, tiny]
        height_ratios = [tiny*aspect_ratio, 1.0 - tiny*aspect_ratio]
        space = 0.05
        hspace = space*aspect_ratio
        wspace = space
        figsize= np.array([aspect_ratio, 1.0])*8

        if self.squishable.binned:
            wavelength = self.squishable.binned_spectral['wavelength']
        else:
            wavelength = self.squishable.spectral['wavelength']
        times = self.squishable.temporal['bjd']
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
        star = self.squishable.target
        self.speak('calculating [{}] for [{}], corrected by the mega-calibrator'.format(key, star))
        target = self.cubes[key][star]
        # KLUDGE!!!
        comparison = self.cubes[key][self.squishable.comparisons[0]]
        z = target/comparison
        oned = np.nanmedian(z, 0)
        z = z/oned[np.newaxis,:]
        return z

    def wavelengthed(self, key='raw_counts', star=None):
        if star == None:
            star = self.squishable.target

        # normalize along the wavelength axis
        self.speak('calculating [{}] for [{}], normalized by its median spectrum'.format(key, star))
        z = self.cubes[key][star]
        oned = np.nanmedian(z, 0)
        return z/oned[np.newaxis,:]


    @property
    def cubes(self):
        if self.squishable.binned:
            return self.squishable.binned_cube
        else:
            return self.squishable.cubes

    def overlayQuality(self, color='tomato'):
        '''
        Plot the cosmic rays on top of the plot.
        '''
        cmap = craftroom.cmaps.one2another(color, color, alphabottom=1.0, alphatop=0.0)
        a = self.loupe.ax['2d']

        try:
            ok = self.ok
        except AttributeError:
            ok = self.cubes['ok'][self.squishable.target]

        self.overlay = a.imshow(ok,
                                    cmap=cmap,
                                    extent=self.loupe.extent,
                                    interpolation='nearest',
                                    zorder=5,
                                    origin='lower',
                                    aspect=self.loupe.ax['2d'].get_aspect()
        )

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
