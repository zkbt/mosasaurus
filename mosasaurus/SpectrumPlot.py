from imports import *
import zachocolor.cmaps

class Color(Talker):
    def __init__(self, range=[400,1050]):
        Talker.__init__(self)
        self.range = range

class Eye(Color):
    def __init__(self, **kwargs):
        Color.__init__(self, **kwargs)

    def color(self, nm):
        return zachopy.color.nm2rgb(nm)

class Gradient(Color):
    def __init__(self, bottom='indigo', top='orange', **kwargs):
        Color.__init__(self, **kwargs)
        self.cmap = zachocolor.cmaps.one2another(bottom, top)
        self.normalizer = plt.Normalize(*self.range)

    def color(self, nm):
        return self.cmap(self.normalizer(np.mean(nm)))



class SpectrumPlot(Talker):
    def __init__(self, xlim=[400,1050]):
        Talker.__init__(self)

        self.toplot = ['spectrum', 'lightcurves', 'depth']#, 'rs_over_a', 'b']

        # populate the plots
        self.panels = {}
        self.panels['spectrum'] = StellarSpectrum(self)
        self.panels['lightcurves'] = LightCurves(self)
        self.panels['depth'] = TransitDepths(self)
        #self.panels['u1'] = FittedParameter(self, key='u1', ylabel='Limb Darkening u1')
        #self.panels['u2'] = FittedParameter(self, key='u2', ylabel='Limb Darkening u2')
        #self.panels['rs_over_a'] = FittedParameter(self, key='rs_over_a', ylabel='Rs/a')
        #self.panels['b'] = FittedParameter(self, key='b', ylabel='b', ylim=[0,1])

        self.xlim = xlim
        self.color=Eye(range=self.xlim).color

    @property
    def npanels(self):
        return len(self.toplot)

    def setup(self, spectrum):

        # set up the gridspec
        plt.figure(figsize=(10,5), dpi=100)
        gs = plt.matplotlib.gridspec.GridSpec(self.npanels, 1,
                        hspace=0.05, wspace=0,
                        left=0.2, bottom=0.15,
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
            ax.set_xlim(*self.xlim)

            if i == (self.npanels-1):
                # set the xlimits
                try:
                    assert(self.xlim[0] is not None)
                except:
                    self.xlim = [(spectrum.bins[0].left - spectrum.binsize/2)/spectrum.unit,
                            (spectrum.bins[-1].right +  spectrum.binsize/2)/spectrum.unit]
                ax.set_xlim(*self.xlim)

                # set the xlabel
                ax.set_xlabel('Wavelength ({0})'.format(spectrum.unitstring))
            else:

                # hide the other xticklabels
                plt.setp(ax.get_xticklabels(), visible=False)
            self.ax[key] = ax

        self.hasbeensetup = True

    def plot(self, spectrum, marker='o', alpha=1.0):

        self.marker=marker
        self.alpha=alpha
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

    def preface(self):
        pass#self.speak('plotting {0} for spectrum {1}'.format(self, self.parent.spectrum))

    def __repr__(self):
        return '<{0}>'.format(self.__class__.__name__)

class StellarSpectrum(Plottable):
    height = 0.5
    ylabel = 'Raw Detected Flux\n(photons/nm/exposure)'

    def plot(self, ax):
        self.preface()

        spectrum = self.parent.spectrum
        wavelength, stellarSpectrum = np.load(spectrum.obs.extractionDirectory + 'medianSpectrum.npy')
        ax.plot(wavelength/spectrum.unit, stellarSpectrum*spectrum.unit, linewidth=3, alpha=0.5, color='black')


class LightCurves(Plottable):
    height = 0.5
    ylabel = 'Time from Mid-Transit\n(hours)'

    def plot(self, ax):
        self.preface()



        spectrum = self.parent.spectrum
        wavelengths = spectrum.wavelengths/spectrum.unit
        depths = spectrum.fitted['k']**2

        # function to normalize lightcurves onto a rotated wavelength scale
        def normalize(flux):
            ratio = 0.8
            one = (flux-1.0)/np.mean(depths)
            return spectrum.binsize/spectrum.unit*ratio*(one+0.5)

        for i in range(len(wavelengths)):
            b = spectrum.bins[i]
            ok = b.tlc.bad == False
            x = normalize(b.tlc.corrected()[ok]) + wavelengths[i]
            y = b.tlc.timefrommidtransit()[ok]

            ax.plot(x, y,
                    marker=self.parent.marker, markeredgewidth=0, linewidth=0, alpha=0.25,
                    color=self.parent.color([b.left/b.unit, b.right/b.unit]))





class FittedParameter(Plottable):
    height = 1.0

    def __init__(self, parent, key=None, ylabel=None, height=1.0, ylim=[None,None]):
        Plottable.__init__(self, parent)
        self.key = key
        self.height = height
        self.ylabel = ylabel
        self.ylim = ylim

    def plot(self, ax):
        self.preface()

        spectrum = self.parent.spectrum
        key = self.key
        parameter = spectrum.fitted[key]
        uncertainty = spectrum.uncertainty[key]
        wavelengths = spectrum.wavelengths/spectrum.unit

        width = 3
        for i in range(len(spectrum.wavelengths)):
            b = spectrum.bins[i]
            ax.errorbar(wavelengths[i], parameter[i], uncertainty[i],
                marker=self.parent.marker, markersize=10, alpha=self.parent.alpha,
                color=self.parent.color([b.left/b.unit, b.right/b.unit]),
                linewidth=width, elinewidth=width, capthick=width, capsize=5, markeredgewidth=0)
        ax.set_ylim(*self.ylim)

class TransitDepths(FittedParameter):
    height = 1.0
    ylabel = 'Transit Depth (%)'
    def __init__(self, parent, ylabel=ylabel, height=height):
        FittedParameter.__init__(self, parent, key='depth', ylabel=ylabel, height=height)

    def plot(self, ax):
        spectrum = self.parent.spectrum
        spectrum.fitted['depth'] = 100*spectrum.fitted['k']**2
        spectrum.uncertainty['depth'] = 100*2*spectrum.uncertainty['k']*spectrum.fitted['k']
        FittedParameter.plot(self, ax)

class FloatingGeometry(SpectrumPlot):
    def __init__(self, **kwargs):
        SpectrumPlot.__init__(self, **kwargs)

        self.toplot = ['spectrum', 'lightcurves', 'depth', 'rs_over_a', 'b']

        # populate the plots
        self.panels = {}
        self.panels['spectrum'] = StellarSpectrum(self)
        self.panels['lightcurves'] = LightCurves(self)
        self.panels['depth'] = TransitDepths(self)
        #self.panels['u1'] = FittedParameter(self, key='u1', ylabel='Limb Darkening u1')
        #self.panels['u2'] = FittedParameter(self, key='u2', ylabel='Limb Darkening u2')
        self.panels['rs_over_a'] = FittedParameter(self, key='rs_over_a', ylabel='Rs/a')
        self.panels['b'] = FittedParameter(self, key='b', ylabel='b', ylim=[0,1])

class FixedGeometry(SpectrumPlot):
    def __init__(self, **kwargs):
        SpectrumPlot.__init__(self, **kwargs)

        self.toplot = ['spectrum', 'lightcurves', 'depth', 'u1', 'u2']

        # populate the plots
        self.panels = {}
        self.panels['spectrum'] = StellarSpectrum(self)
        self.panels['lightcurves'] = LightCurves(self)
        self.panels['depth'] = TransitDepths(self)
        self.panels['u1'] = FittedParameter(self, key='u1', ylabel='Limb Darkening u1')
        self.panels['u2'] = FittedParameter(self, key='u2', ylabel='Limb Darkening u2')

class Combined(SpectrumPlot):
    def __init__(self, **kwargs):
        SpectrumPlot.__init__(self, **kwargs)

        self.toplot = ['depth', 'u1', 'u2']

        # populate the plots
        self.panels = {}
        self.panels['depth'] = TransitDepths(self)
        self.panels['u1'] = FittedParameter(self, key='u1', ylabel='Limb Darkening u1')
        self.panels['u2'] = FittedParameter(self, key='u2', ylabel='Limb Darkening u2')

class JustDepth(SpectrumPlot):
    def __init__(self, **kwargs):
        SpectrumPlot.__init__(self, **kwargs)

        self.toplot = ['depth']

        # populate the plots
        self.panels = {}
        self.panels['depth'] = TransitDepths(self)
