import TransmissionSpectrum
t = TransmissionSpectrum.TransmissionSpectrum('../wasp94_140809.obs')
b = t.bins[0]
b.fit(*t.psi(), slow=True, plot=True, interactive=True, nburnin=50, ninference=50)
