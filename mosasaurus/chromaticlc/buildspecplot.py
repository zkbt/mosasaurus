from imports import *
import CombinedTransmissionSpectrum
c = CombinedTransmissionSpectrum.CombinedTransmissionSpectrum(binsize=50)
c.load()
c.combine()

markers=['^', 's', 'D']
nights=['1 Aug 2014', '5 Aug 2014', '9 Aug 2014']
import SpectrumPlot
reload(SpectrumPlot)
s = SpectrumPlot.JustDepth(xlim=[440,1000])
filename = ''
for i in range(len(c.spectra)):
    t = c.spectra[i]
    s.plot(t, alpha=0.2, marker=markers[i])
    s.ax['depth'].set_ylim(1,1.3)
    s.ax['depth'].set_title('WASP94Ab Transmission Spectrum')
    filename += "{0}".format(i)

    plt.savefig(filename + '.pdf')
s.plot(c)
i += 1
filename += "{0}".format(i)
plt.savefig(filename + '.pdf')
