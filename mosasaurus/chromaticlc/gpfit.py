# for testing, simply fit one bin
whichbin = 0
# how wide should the bins be?
size = 100
binsize = 5.0
import zachopy.utils, transit, matplotlib.pyplot as plt, matplotlib

# the initial conditions are stored in their own script
import initialConditions

# create the combined transmission spectrum object
from mosasaurus.CombinedTransmissionSpectrum import CombinedTransmissionSpectrum
c = CombinedTransmissionSpectrum(
    files=['wasp94_140801.obs', 'wasp94_140805.obs', 'wasp94_140809.obs'],
    binsize=size)

# run a pure LM fit, using the old LM fit
kw = dict(  remake=False, verbose=False,
            wobbly=False, floatLD=False,
            mcmc=False)

for w in c.wavelengths:
    c.fitMonochromatic(w, initialConditions, plot=False, **kw)
    lm = c.synthesizer
    for tlc in lm.tlcs:
        tlc.binto(5.0/60.0/24.0)
    lm.trainGP(plot=True)

    mcmc = transit.MCMC(tlcs=lm.tlcs,
                directory=c.synthesizer.directory.replace('/lm/', '/mcmc/'),
                likelihoodtype='red_gp')
    for thing in ['period', 't0',  'b', 'rsovera', 'semiamplitude']:
        mcmc.tieAcrossEpochs(thing)
        mcmc.tieAcrossTelescopes(thing)

    for thing in ['k', 'u1', 'u2', 'gplna', 'gplntau']:
        mcmc.tieAcrossEpochs(thing)
    mcmc.fit(fromcovariance=False, nwalkers=200, updates=3, plot=True, **kw)


# print the fitted parameters
c.synthesizer.pdf.printParameters()
