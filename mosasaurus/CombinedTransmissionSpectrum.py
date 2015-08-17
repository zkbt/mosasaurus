from imports import *
from TransmissionSpectrum import TransmissionSpectrum
from WavelengthBin import WavelengthBin

class CombinedObs(Talker):

    def __init__(self, name):
        self.name = name
        self.night = 'combined'

class CombinedTransmissionSpectrum(TransmissionSpectrum):
    def __init__(self, method='lm', label='fixedGeometry', files=['../wasp94_140801.obs', '../wasp94_140805.obs', '../wasp94_140809.obs'], binsize=100):
        Talker.__init__(self)


        # keep track of the files
        self.files = files

        # populate list of transmission spectra
        self.spectra = []
        self.speak('creating a combined transmission spectrum from files:')
        for f in files:
            self.speak('  {0}'.format(f))
            spectrum = TransmissionSpectrum(f, binsize=binsize)
            self.spectra.append(spectrum)

        self.obs = CombinedObs(self.spectra[0].obs.name)

        self.setBinsize(self.spectra[0].binsize)

        self.mask = 'defaultMask'
        self.label = label

        binoptions = []
        for s in self.spectra:
            for b in s.bins:
                binoptions.append(b.__repr__().split('|')[-1].split('>')[0])
                print binoptions

        bins, wavelengths = [], []
        for b in np.unique(binoptions):
            left, right = np.array(b.split('nm')[0].split('to')).astype(np.int)*10
            bins.append(WavelengthBin(self, left=left, right=right))
            wavelengths.append(bins[-1].wavelength)

        s = np.argsort(wavelengths)
        self.bins = np.array(bins)[s]
        self.wavelengths = np.array(wavelengths)[s]

    def load(self, method='lm'):
        self.speak('loading all spectra')
        for s in self.spectra:
            s.load(method)

    def combine(self):

        self.fitted, self.uncertainty, self.chisq, self.dof = {}, {}, {}, {}

        # loop over wavelengths

        for i in range(len(self.wavelengths)):
            w = self.wavelengths[i]
            self.speak('combining all available spectra for {0}'.format(w))
            f2combine, u2combine = {},{}

            for s in self.spectra:
                ok = (s.wavelengths == w).nonzero()
                self.speak('wavelength {0} is bin {1}'.format(w, ok))
                assert(len(ok) == 1)
                for k in s.fitted.keys():
                    try:
                        f2combine[k], u2combine[k]
                    except KeyError:
                        f2combine[k], u2combine[k] = [], []
                    f2combine[k].append(s.fitted[k][ok][0])
                    u2combine[k].append(s.uncertainty[k][ok][0])

            for k in f2combine.keys():
                f2combine[k] = np.array(f2combine[k])
                u2combine[k] = np.array(u2combine[k])
                try:
                    self.fitted[k], self.uncertainty[k]
                except:
                    self.fitted[k] = np.empty(len(self.wavelengths))
                    self.uncertainty[k] = np.empty(len(self.wavelengths))
                    self.chisq[k] = np.empty(len(self.wavelengths))
                    self.dof[k] = np.empty(len(self.wavelengths))

                self.fitted[k][i] = np.sum(f2combine[k]/u2combine[k]**2)/np.sum(1.0/u2combine[k]**2)
                self.chisq[k][i] = np.sum((f2combine[k] - self.fitted[k][i])**2/u2combine[k]**2)
                self.dof[k][i] = max(len(f2combine[k]) - 1.0, 1)
                rescaling = max(np.sqrt(self.chisq[k][i]/self.dof[k][i]), 1)

                self.uncertainty[k][i] = 1.0/np.sqrt(np.sum(1.0/u2combine[k]**2))*rescaling

            print f2combine
            print u2combine
            print

                ######## PICK UP HERE!

        pass

    def table(self):
        form = '{0:>20}{1:>20}{2:>20}{3:>25}'
        print form.format('left', 'right', 'rp_over_rs', 'rp_over_rs_error')

        for i in range(len(self.bins)):
            b = self.bins[i]
            if i == 0:
                print form.format(b.unitstring, b.unitstring, 'unity', 'unity')
            print form.format(b.left/b.unit, b.right/b.unit, self.fitted['k'][i], self.uncertainty['k'][i])

    @property
    def mask(self):
        "Indicate the mask being using for all spectra."
        return self._mask

    @mask.setter
    def mask(self, value):
        "Specify the mask to use (for all spectra)."

        self.speak('setting mask for all spectra to {0}'.format(value))
        # update the maskname hidden value
        self._mask = value

        # loop over the transmission spectra
        for spectrum in self.spectra:
            self.speak('setting mask for {spectrum} to {mask}'.format(spectrum=spectrum, mask=value))
            # use the transmission spectrum's applyMask method to load the appropriate mask and apply it to the light curves
            spectrum.applyMask(value)

    @property
    def label(self):
        "Indicate the fit being used for all spectra."
        return self._label

    @label.setter
    def label(self, value):
        "Specify the mask to use (for all spectra)."

        self.speak('setting label for all spectra to {0}'.format(value))
        self._label = value
        for spectrum in self.spectra:
            self.speak('setting label for {spectrum} to {label}'.format(spectrum=spectrum, label=value))
            spectrum.label = value

    def __repr__(self):
        """How should this object be represented (e.g. when shown as an element in a list)"""
        s = "<CombinedTransmissionSpectrum "
        for t in self.spectra:
            s += '{0}|'.format( repr(t))
        s = s.strip('|') +  '>'
        return s
