from CombinedTransmissionSpectrum import *

class CombinedEmissionSpectrum(CombinedTransmissionSpectrum):
    def __init__(self, *args, **kwargs):
        self.phaseofinterest = 0.5
        CombinedTransmissionSpectrum.__init__(self, *args, **kwargs)

    def fitMonochromatic(self,
                            w, # wavelength dictionary key to fit
                            initialConditions=None, # required init. conditions
                            likelihoodtype='white',
                            mcmc=False,
                            remake=False,
                            plot=True,
                            verbose=False,
                            label=None,
                            ploteverything=True,
                            **kw):

        ''' fit a wavelength bin, across all observations,
            given some initial conditions '''

        likelihoodtypes = ['white', 'red_beta', 'red_gp']
        try:
            assert(likelihoodtype in likelihoodtypes)
        except:
            raise ValueError(
                "likelihoodtype is not in {0}".format(likelihoodtypes))

        # create an empty list of tlcs
        tlcs = []

        # keep track of labels for this fit
        if label is None:
            label='eclipse'
        self.label = label

        # loop over the tlcs, and create a model for each
        for i, orig in enumerate(self.archiveoftlcs[w]):


            # define the input objects
            planet = transit.Planet(**initialConditions.planetkw)

            # create a synthesizer director for this bin
            synthesizerdirectory = '{0}{1}/'.format(
                                        self.fitdirectory,
                                        self.bins[self.w2bin(w)][0].identifier
                                        )
            zachopy.utils.mkdir(synthesizerdirectory)

            # assign an epoch to the TLC
            tlc = orig.splitIntoEpochs(planet, newdirectory=synthesizerdirectory, phaseofinterest=self.phaseofinterest)[0]

            # store that epoch in the archive of light curves
            self.archiveoftlcs[w][i] = tlc

            # also, add it to the list to be included in this fit
            tlcs.append(tlc)

            # turn verbosity on or off for the input TLCs
            if verbose:
                tlc.pithy=False
            else:
                tlc.pithy=True

            # float the planetary parameters (including geometry, if desired)
            planet.J.float(limits=[-0.05, 0.05])#(limits=[-0.1, 0.1])
            planet.dt.float(limits=[-0.07, 0.07])

            # float the GP hyperparameters, or simply the normal linear basis functions
            #if likelihoodtype == 'gp':
            instrument = transit.Instrument(tlc=tlc, gplna=-5, gplntau=-5, **initialConditions.instrumentkw)
                #instrument.gplna.float(-10,[-20,0])
                #instrument.gplntau.float(-5, [-10,0])


            # a constant baseline
            instrument.C.float(value=1.0,limits=[0.9, 1.1])
            # instrument rotator angle (seems to matter!)
            instrument.rotatore_tothe1.float(value=0.002, limits=[-0.005, 0.005])


            # width of the whole spectrum, and the width in the wavelength range
            instrument.width_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])
            instrument.dwidth_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])

            # cross-dispersion centroid of the whole spectrum, and in the wavelength range
            instrument.centroid_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])
            instrument.dcentroid_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])

            # applied wavelength offset
            try:
                instrument.shift_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])
            except AttributeError:
                self.speak("couldn't find shifts!")

            # the sky brightness in
            instrument.sky_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])
            instrument.peak_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])

            star = transit.Star()
            # create the transit model
            tm = transit.TM(planet, star, instrument, directory=tlc.directory)
            tm.linkLightCurve(tlc)

        if mcmc:
            self.synthesizer = transit.MCMC(tlcs=tlcs, directory=synthesizerdirectory, likelihoodtype=likelihoodtype)
        else:
            self.synthesizer = transit.LM(tlcs=tlcs, directory=synthesizerdirectory)

        for thing in []:
            self.synthesizer.tieAcrossEpochs(thing)
            self.synthesizer.tieAcrossTelescopes(thing)

        for thing in ['J', 'gplna', 'gplntau']:
            self.synthesizer.tieAcrossEpochs(thing)

        self.synthesizer.speak('the starting parameters are')
        self.synthesizer.printParameters()

        try:
            self.archiveofpdfs
        except AttributeError:
            self.archiveofpdfs = {}

        self.synthesizer.fit(remake=remake, fromcovariance=False, **kw)
        if mcmc:
            pass
        else:
            self.synthesizer.pdf.calculateUncertainties(style='covariance')
        self.archiveofpdfs[w] = self.synthesizer.pdf

        if plot:
            transit.MultiplexPlot(self.synthesizer.tlcs, transit.DiagnosticsPlots,
                                wobbly=True, everything=ploteverything,
                                figsize=(30,10), dpi=72, ylim=[0.995,1.005],
                                binsize=15/60./24.0, phaseofinterest=0.5)
            plt.savefig(self.synthesizer.directory + 'lightcurves.pdf')
