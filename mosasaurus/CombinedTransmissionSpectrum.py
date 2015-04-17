from imports import *
from TransmissionSpectrum import TransmissionSpectrum


class CombinedTransmissionSpectrum(TransmissionSpectrum):
    def __init__(self, files=['../wasp94_140801.obs', '../wasp94_140805.obs', '../wasp94_140809.obs']):
        Talker.__init__(self)


        # keep track of the files
        self.files = files

        # populate list of transmission spectra
        self.spectra = []
        self.speak('creating a combined transmission spectrum from files:')
        for f in files:
            self.speak('  {0}'.format(f))
            spectrum = TransmissionSpectrum(f)
            spectrum.loadLCs()
            self.spectra.append(spectrum)

        self.mask = 'defaultMask'
        self.label = 'fixedGeometry'

    @property
    def mask(self):
        "Indicate the mask being using for all spectra."
        return self._mask

    @mask.setter
    def mask(self, value):
        "Specify the mask to use (for all spectra)."

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
