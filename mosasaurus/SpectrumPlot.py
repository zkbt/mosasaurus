from imports import *

class SpectrumPlot(Talker):
    def __init__(self):
        Talker.__init__(self)

        self.toplot = ['spectrum', 'depth']

        # populate the plots
        self.panels = {}
        self.panels['spectrum'] = StellarSpectrum(self)
        #self.panels['lightcurves'] = LightCurves(self)
        self.panels['depth'] = TransitDepths(self)

    @property
    def npanels(self):
        return len(self.toplot)

    def setup(self, spectrum, xlim=[400,1050]):

        # set up the gridspec
        gs = plt.matplotlib.gridspec.GridSpec(self.npanels, 1,
                        hspace=0.1, wspace=0,
                        left=0.2,
                        height_ratios=[self.panels[k].height for k in self.toplot])

        # create the axes for all the plots
        sharex=None
        self.ax = {}
        for i in range(self.npanels):
            # pull out keys in the order specified in toplot
            key = self.toplot[i]

            # create the axes object, and make sure its xaxis is shared
            ax = plt.subplot(gs[i], sharex=sharex)
            sharex = ax

            p = self.panels[key]
            ax.set_ylabel(p.ylabel)
            ax.set_xlim(*xlim)

            if i == (self.npanels-1):
                # set the xlimits
                if xlim is None:
                    xlim = [(spectrum.bins[0].left - spectrum.binsize/2)/spectrum.unit,
                            (spectrum.bins[-1].right +  spectrum.binsize/2)/spectrum.unit]
                ax.set_xlim(*xlim)

                # set the xlabel
                ax.set_xlabel('Wavelength ({0})'.format(spectrum.unitstring))
            else:

                # hide the other xticklabels
                plt.setp(ax.get_xticklabels(), visible=False)
            self.ax[key] = ax

        self.hasbeensetup = True

    def plot(self, spectrum):

        try:
            assert(self.hasbeensetup)
        except:
            self.setup(spectrum)

        # point the plot at this particular transmission spectrum
        self.spectrum = spectrum

        # loop through the keys, and plot them
        for k in self.toplot:
            panel = self.panels[k]
            self.speak('adding {self.spectrum} to {panel}'.format(**locals()))
            panel.plot(self.ax[k])

class Plottable(Talker):
    def __init__(self, parent):
        Talker.__init__(self)
        self.parent = parent


class StellarSpectrum(Plottable):
    height = 0.5
    ylabel = 'Raw Detected Flux\n(photons/nm/exposure)'

    def plot(self, ax):
        # plot the spectrum
        if True:
            spectrum = self.parent.spectrum
            wavelength, stellarSpectrum = np.load(spectrum.obs.extractionDirectory + 'medianSpectrum.npy')
            ax.plot(wavelength/spectrum.unit, stellarSpectrum*spectrum.unit, linewidth=3, alpha=0.5, color='black')


class LightCurves(Plottable):
    height = 0.5
    ylabel = 'Time from Mid-Transit\n(hours)'

    def plot(self, ax):
        pass

class TransitDepths(Plottable):
    height = 1.0
    ylabel = 'Transit Depth (%)'

    def plot(self, ax):
        pass

        # function to normalize lightcurves onto a rotated wavelength scale
        def normalize(flux):
            ratio = 0.8
            one = (flux-1.0)/np.mean(self.fitted['k'])**2
            return self.binsize/self.unit*ratio*(one+0.5)



def plot(self):
  self.setupSuperPlot()
  colors = []



  for i in np.arange(len(self.bins)):
    # select this bin
    bin = self.bins[i]

    # plot the model for this bin
    #time, planetmodel, instrumentmodel = bin.tlc.TM.smooth_model()
    #kw = {'color':zachopy.color.nm2rgb([bin.left/self.unit, bin.right/self.unit], intensity=0.25), 'linewidth':3, 'alpha':0.5}
    #self.ax_lc.plot(normalize(planetmodel) + bin.wavelength/self.unit, bin.tlc.TM.planet.timefrommidtransit(time), **kw)

    # plot the (instrument corrected) datapoints
    kw = {'marker':'.', 'color':zachopy.color.nm2rgb([bin.left/self.unit, bin.right/self.unit]), 'alpha':0.25, 'linewidth':0, 'marker':'o'}
    self.ax_lc.plot(normalize(bin.tlc.corrected()[bin.tlc.bad == False]) + bin.wavelength/self.unit, bin.tlc.timefrommidtransit()[bin.tlc.bad == False], **kw)
    colors.append(kw['color'])

    print bin
    print bin.tlc.TM.planet

    width = 3

    self.ax_ts.errorbar(self.wavelengths[i]/self.unit, self.fitted['k'][i], self.uncertainty['k'][i], marker='o', color=zachopy.color.nm2rgb([bin.left/self.unit, bin.right/self.unit]), markersize=10, linewidth=width, elinewidth=width, capsize=5, capthick=width)
