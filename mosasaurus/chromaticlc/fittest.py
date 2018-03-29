import TransmissionSpectrum
t = TransmissionSpectrum.TransmissionSpectrum('../wasp94_140809.obs', binsize=50)
t.applyMask('defaultMask')
t.fitBins(slow=False, plot=True, label='fixedGeometry', remake=True)


import TransmissionSpectrum
t = TransmissionSpectrum.TransmissionSpectrum('../wasp94_140809.obs')
t.applyMask('defaultMask')
t.fitBins(slow=False, plot=True, label='fixedGeometry', remake=True)



import SpectrumPlot
s = SpectrumPlot.FloatingGeometry()
nights = [1,5,9]
markers = {1:'o', 5:'s', 9:'^'}
for night in nights:
    file = '../wasp94_14080{0}.obs'.format(night)
    t = TransmissionSpectrum.TransmissionSpectrum(file)
    t.label='floatingGeometry'
    t.applyMask('defaultMask')
    t.load('lm')
    s.plot(t, marker=markers[night])


import SpectrumPlot
s = SpectrumPlot.FixedGeometry()
nights = [1,5,9]
markers = {1:'o', 5:'s', 9:'^'}
for night in nights:
    file = '../wasp94_14080{0}.obs'.format(night)
    t = TransmissionSpectrum.TransmissionSpectrum(file)
    t.label='fixedGeometry'
    t.applyMask('defaultMask')
    t.load('lm')
    s.plot(t, marker=markers[night])



#t.fitBins(slow=True, plot=True, nburnin=500, ninference=50, label='floatingGeometry')
#b = t.bins[0]
#b.fit(slow=False, plot=True, remake=True)

#b.fit(slow=True, plot=True, interactive=True, nburnin=50, ninference=50)
