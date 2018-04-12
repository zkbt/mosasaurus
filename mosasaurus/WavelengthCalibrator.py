'''Wavelength Calibrator defines the wavelength calibration, and '''

from .imports import *
from numpy.polynomial import Legendre
colors = dict(He='darkorange', Ne='red', Ar='deepskyblue')
shortcuts = {'h':'He', 'n':'Ne', 'a':'Ar'}

class WavelengthCalibrator(Talker):
    def __init__(self,  aperture,
                        elements=['He', 'Ne','Ar'],
                        polynomialdegree=3,
                        matchdistance=100):
        Talker.__init__(self)

        # keep track of the aperture this belongs to
        self.aperture = aperture

        # set up the some defaults
        self.elements = elements
        self.polynomialdegree = polynomialdegree
        self.matchdistance = matchdistance

        # either load a previous wavecal, or create a new one
        self.populate()

    @property
    def wavelengthprefix(self):
        '''Create the start of a filename to store wavelength calibration information.'''
        return os.path.join(self.aperture.directory, '{0}_'.format(self.aperture.name))

    @property
    def calibrationfilename(self):
        '''Filename to store the wavelength calibration.'''
        return self.wavelengthprefix + 'wavelengthcalibration.npy'

    @property
    def waveidfilename(self):
        '''Filename to store the identified wavelength matches for this solution.'''
        return self.wavelengthprefix + 'waveids.txt'

    def loadWavelengthIdentifications(self, restart=False):
        ''' Try loading a custom stored wavelength ID file
            from the aperture's directory, and if that
            doesn't work, load the default one for this grism.'''

        try:
            # try to load a custom wavelength id file
            assert(restart == False)

            self.speak('checking for a custom wavelength-to-pixel file')
            self.waveids = astropy.io.ascii.read(self.waveidfilename)

            # keep track of which waveid file is used
            self.whichwaveid = 'Aperture-Specific ({0})'.format(
                self.waveidfilename.split('/')[-1])

        except (IOError,AssertionError):
            self.speak('no custom wavelength-to-pixel files found')

            # load the default for this grism, as set in obs. file
            d = astropy.io.ascii.read(self.aperture.instrument.wavelength2pixelsFile)
            self.rawwaveids = d[['pixel', 'wavelength', 'name']]
            self.rawwaveids['pixel'] /= self.aperture.instrument.binning

            # use a cross-corrlation to find the rough offset
            #  (the function call will define waveids)
            # pull out the peaks from extracted arc spectra
            self.findPeaks()

            # perform a first rough alignment, to make pixels
            self.findRoughShift()


            # keep track of whcih waveid file is used
            self.whichwaveid = 'Default ({0})'.format(
                self.aperture.instrument.wavelength2pixelsFile.split('/')[-1])


        self.speak('loaded {0}:'.format(self.whichwaveid))
        self.findKnownWavelengths()
        self.guessMatches()

    def saveWavelengthIdentification(self):
        '''store the wavelength-to-pixel identifications'''

        self.waveids.write( self.waveidfilename,
                            format='ascii.fixed_width',
                            delimiter='|',
                            bookend=False)
        self.speak('saved wavelength-to-pixel matches to '+self.waveidfilename)

    def save(self):

        self.speak('saving all aspects of the wavelength calibration')

        # the wavelength identifications
        self.saveWavelengthIdentification()

        # the actual matches used (with outliers, etc...)
        self.saveMatches()

        # save the coefficients (and domain) of the calibration
        self.saveCalibration()

        # save the figure
        self.figcal.savefig(self.wavelengthprefix + 'calibration.pdf')

    def loadCalibration(self):

        coef, domain = np.load(self.calibrationfilename)
        self.pixelstowavelengths = Legendre(coef, domain)
        self.polynomialdegree = self.pixelstowavelengths.degree()
        self.speak("loaded wavelength calibration"
                "from {0}".format(self.calibrationfilename))

    def saveCalibration(self):
        np.save(self.calibrationfilename, (self.pixelstowavelengths.coef, self.pixelstowavelengths.domain))
        self.speak("saved wavelength calibration coefficients to {0}".format(self.calibrationfilename))

    def populate(self, restart=False):

        # populate the wavelength identifications
        self.loadWavelengthIdentifications(restart=restart)
        try:
            # populate the wavelength calibration polynomial
            assert(restart==False)
            self.loadCalibration()
            self.justloaded = True
            self.loadMatches()
            #self.plotWavelengthFit(interactive=False)
            #unhappy = ('n' in self.input('Are you happy with the wavelength calibration? [Y,n]').lower())
            #assert(unhappy == False)
        except (IOError, AssertionError):
            self.justloaded = False
            self.create()


    def findRoughShift(self, blob=2.0):
        '''using a list of pixels matched to wavelengths, find the rough
                offset of the this arc relative to the standard slits/setting'''

        self.speak("cross correlating arcs with known wavelengths")

        # create a plot showing how well the lines match
        '''
        figure_waverough = plt.figure(  'wavelength rough offset',
                                        figsize=(6,4), dpi=100)
        gs = plt.matplotlib.gridspec.GridSpec(  3,2,
                                                bottom=0.15, top=0.85,
                                                hspace=0.1,wspace=0,
                                                width_ratios=[1, .5])

        # make the axes for checking the rough alignment
        self.ax_waverough, self.ax_wavecor = {}, {}
        sharer, sharec = None, None
        for i, e in enumerate(self.elements):
            self.ax_waverough[e] = plt.subplot(gs[i,0],sharex=sharer)
            self.ax_wavecor[e] = plt.subplot(gs[i,1], sharex=sharec)
            sharer, sharec = self.ax_waverough[e], self.ax_wavecor[e]

        # calculate correlation functions
        self.corre = {}
        for count, element in enumerate(self.elements):
            # pull out the peaks
            xPeak = [p['w'] for p in self.peaks[element]]
            yPeak = [p['intensity'] for p in self.peaks[element]]

            # create fake spectra using the line positions (reference + new)
            x = np.arange(-self.aperture.instrument.ysize,self.aperture.instrument.ysize)

            myPeaks, theirPeaks = np.zeros(len(x)), np.zeros(len(x))
            # create fake spectrum of their peaks
            for i in range(len(self.rawwaveids)):
                if element in self.rawwaveids['name'][i]:
                    center = self.rawwaveids['pixel'][i]
                    theirPeaks += np.exp(-0.5*((x-center)/blob)**2)

            # create fake spectrum of my peaks
            for i in range(len(xPeak)):
                center = xPeak[i]
                myPeaks += np.exp(-0.5*((x-center)/blob)**2)*np.log(yPeak[i])

            # calculate the correlation function for this element
            self.corre[element] = np.correlate(myPeaks, theirPeaks, 'full')

            # plot the rough shift and identifications
            self.ax_waverough[element].plot(x, myPeaks/myPeaks.max(),
                                            label='extracted', alpha=0.5,
                                            color=colors[element])
            self.ax_waverough[element].set_ylim(0, 2)


            # plot the correlation functions
            normcor = self.corre[element]
            assert(np.isfinite(self.corre[element]).any())

            normcor /= np.nanmax(self.corre[element])
            self.ax_wavecor[element].plot(normcor, label=element, alpha=0.5,
                                            color=colors[element])
            # tidy up the plots
            for a in [self.ax_wavecor, self.ax_waverough]:
                plt.setp(a[element].get_xticklabels(), visible=False)
                plt.setp(a[element].get_yticklabels(), visible=False)

            # multiply the correlation functions together
            assert(np.isfinite(self.corre[element]).any())
            if count == 0:
                self.corre['combined'] = np.ones_like(self.corre[element])
            self.corre['combined'] *= self.corre[element]
        estimatedpeakoffset = np.argmax(self.corre['combined']) - len(x)
        '''

        # define the new, shifted, waveids array
        self.waveids = copy.deepcopy(self.rawwaveids)
        self.waveids['pixel'] += self.aperture.obs.instrument.peakoffset

        '''
        # plot the shifted wavelength ids, and combined corfuncs
        for element in self.elements:
            for i in range(len(self.rawwaveids)):
                if element in self.rawwaveids['name'][i]:
                      center = self.waveids['pixel'][i]
                      self.ax_waverough[element].axvline(center, alpha=0.25, color='black')
            # plot the combined correlation function
            normedcombined = self.corre['combined']/np.max(self.corre['combined'])
            self.ax_wavecor[element].plot(normedcombined,
                                label='combined', alpha=0.25, color='black')
            # tidy up the plots
            self.ax_waverough[element].set_ylabel(element)
            ax = self.ax_waverough[self.elements[-1]]
            plt.setp(ax.get_xticklabels(), visible=True)
            ax.set_xlabel('Pixel Position')
            fontsize = 8
            self.ax_wavecor[self.elements[0]].set_title('cross correlation peaks at \n{0} pixels ({1}x{1} binned pixels)'.format(estimatedpeakoffset, self.aperture.instrument.binning), fontsize=fontsize)
            self.ax_waverough[self.elements[0]].set_title(
            'Coarse Wavelength Alignment\nfor ({0:0.1f},{1:0.1f})'.format(
                 self.aperture.x, self.aperture.y),fontsize=fontsize)

            # save the figure
            figure_waverough.savefig(
                os.path.join(self.aperture.directory, 'roughWavelengthAlignment_{0}.pdf'.format(
                                self.aperture.name)))
        '''
    @property
    def peaks(self):
        try:
            return self._peaks
        except AttributeError:
            self.findPeaks()
            return self._peaks

    def findPeaks(self):
        '''identify peaks in the extracted arc spectrum'''

        # extract a spectrum from the master image for each lamp
        self.aperture.arcs = {}
        for element in self.elements:
            self.aperture.arcs[element] = self.aperture.extract(
                                n=element,
                                image=self.aperture.images[element],
                                arc=True)



        # find the peaks in my spectra
        self._peaks = {}

        for count, element in enumerate(self.elements):

            # the pixel spectrum self.aperture.extracted from this arc lamp
            width = np.min(self.aperture.trace.extractionwidths)
            flux = self.aperture.arcs[element][width]['raw_counts']

             # identify my peaks
            xPeak, yPeak, xfiltered, yfiltered = craftroom.oned.peaks(
                                                self.aperture.waxis,
                                                flux,
                                                plot=False,
                                                xsmooth=30,
                                                threshold=100,
                                                edgebuffer=10,
                                                widthguess=1,
                                                maskwidth=3,
                                                returnfiltered=True)
            self.aperture.arcs[element]['filtered'] = yfiltered

            # for some reason, need to trim peaks outside range
            pad = 25
            toright = xPeak > (np.min(self.aperture.waxis) + pad)
            toleft = xPeak < (np.max(self.aperture.waxis) - pad)
            ok = toleft*toright
            xPeak, yPeak = xPeak[ok], yPeak[ok]

            n = len(xPeak)

            # store those peaks
            self._peaks[element] = []
            for i in range(n):
                peak =  {
                        'w':xPeak[i],
                        'intensity':yPeak[i],
                        'handpicked':False,
                        'outlier':False
                        }
                self._peaks[element].append(peak)

    def findKnownWavelengths(self):
        # create a temporary calibration to match reference wavelengths to reference pixels (so we can extrapolate to additional wavelengths not recorded in the dispersion solution file)




        self.knownwavelengths = {}

        # treat the arc lamps separately
        for count, element in enumerate(self.elements):
            # pull out the wavelengths from the complete file
            self.knownwavelengths[element] = []
            for i in range(len(self.waveids)):
                if element in self.waveids['name'][i]:
                    wave = self.waveids['wavelength'][i]
                    pixel = self.waveids['pixel'][i]
                    known = {
                            'element':element,
                            'wavelength':wave,
                            'pixelguess':pixel
                            }
                    self.knownwavelengths[element].append(known)

    @property
    def matchesfilename(self):
        return self.wavelengthprefix + 'wavelengthmatches.npy'

    def saveMatches(self):
        self.speak('saving wavelength dictionaries to {}'.format(self.matchesfilename))
        np.save(self.matchesfilename, (self.matches, self.knownwavelengths))

    def loadMatches(self):
        (self.matches, self.knownwavelengths) = np.load(self.matchesfilename)
        self.speak('loaded wavelength matches from {0}'.format(self.matchesfilename))

    def guessMatches(self):
        self.matches = []

        # do identification with one arc at a time
        for element in self.elements:
            # pull out my peaks and theirs
            myPeaks = np.array([p['w'] for p in self.peaks[element]])
            theirPeaksOnMyPixels = np.array([known['pixelguess']
                                            for known
                                            in self.knownwavelengths[element]])

            thiselement = []
            # loop over my peaks
            for i in range(len(myPeaks)):

                # find my closest peak to theirs
                distance = myPeaks[i] - theirPeaksOnMyPixels
                closest = np.sort(np.nonzero(np.abs(distance) == np.min(np.abs(distance)))[0])[0]

                if distance[closest] < self.matchdistance:
                    # call this a match (no matter how far)
                    match = {
                            'mine':self.peaks[element][i],
                            'theirs':self.knownwavelengths[element][closest],
                            'distance':distance[closest]
                            }
                    thiselement.append(match)

            self.matches.extend(thiselement)
            '''for theirs in theirPeaksOnMyPixels:


                relevant = [m for m in thiselement if
                            m['theirs']['pixelguess'] == theirs]

                if len(relevant) > 0:
                    print("{0} of my peaks match to their {1}".format(len(relevant), theirs))

                    distances = np.abs([m['theirs']['pixelguess'] - m['mine']['w'] for m in relevant])


                    best = distances.argmin()

                    print("the closests one is {0}".format(relevant[best]['mine']['w']))
                    self.matches.append(relevant[best])
            '''



            # add this to the list
            #self.matches.extend(thiselement)

    def create(self, remake=False):
        '''Populate the wavelength calibration for this aperture.'''

        self.speak("populating wavelength calibration")



        # loop through
        self.notconverged = True
        while(self.notconverged):

            # set an initial guess matching wavelengths to my pixles
            #self.guessMatches()

            # do a fit
            if self.justloaded:
                self.justloaded = False
            else:
                # do an initial fit
                self.pixelstowavelengths = Legendre.fit(
                                                    x=self.pixel,
                                                    y=self.wavelength,
                                                    deg=self.polynomialdegree,
                                                    w=self.weights
                                                    )
            # identify outliers
            limit = 1.48*craftroom.oned.mad(self.residuals[self.good])*4
            limit = np.maximum(limit, 1.0)
            # outliers get reset each time (so don't lose at edges)
            outlier = np.abs(self.residuals) > limit
            # keep track of which are outlier
            for i, m in enumerate(self.matches):
                self.matches[i]['mine']['outlier'] = outlier[i]
            self.speak('points beyond {0} are outliers ({1})'.format(limit, i))

            self.pixelstowavelengths = Legendre.fit(
                                                    x=self.pixel,
                                                    y=self.wavelength,
                                                    deg=self.polynomialdegree,
                                                    w=self.weights)

            self.plotWavelengthFit()
            #self.updatew2p()

    @property
    def weights(self):
        return self.good+self.handpicked*10



    @property
    def residuals(self):
        return  self.wavelength - self.pixelstowavelengths(self.pixel)

    @property
    def good(self):
        isntoutlier = np.array([m['mine']['outlier'] == False for m in self.matches])
        ishandpicked = np.array([m['mine']['handpicked'] for m in self.matches])
        return (isntoutlier + ishandpicked) > 0

    @property
    def pixel(self):
        return np.array([m['mine']['w'] for m in self.matches])

    @property
    def pixelelement(self):
        return np.array([m['theirs']['element'] for m in self.matches])

    @property
    def intensity(self):
        return np.array([m['mine']['intensity'] for m in self.matches])

    @property
    def wavelength(self):
        return np.array([m['theirs']['wavelength'] for m in self.matches])

    @property
    def emissioncolor(self):
        return np.array([colors[m['theirs']['element']] for m in self.matches])

    def plotWavelengthFit(self, interactive=True):
        # plot to make sure the wavelength calibration makes sense

        self.speak('{}, {}'.format(self.pixelstowavelengths, self.pixelstowavelengths.domain, self.pixelstowavelengths.window))
        self.figcal = plt.figure('wavelength calibration',
                                        figsize=(15,6), dpi=72)
        self.interactivewave = craftroom.displays.iplot.iplot(4,1,
                height_ratios=[0.1, 0.4, 0.2, .2], hspace=0.1,
                bottom=0.15)

        self.ax_header = self.interactivewave.subplot(0)

        # do the lamp spectra overlap?
        self.ax_walign = self.interactivewave.subplot(1)

        # what is the actual wavelength calibration
        self.ax_wcal = self.interactivewave.subplot(2, sharex=self.ax_walign)

        # what are the residuals from the fit
        self.ax_wres = self.interactivewave.subplot(3, sharex=self.ax_walign)

        for ax in [self.ax_header, self.ax_walign, self.ax_wcal]:
            plt.setp(ax.get_xticklabels(), visible=False)

        self.ax_header.set_title("Wavelength Calib. for Aperture"
            " (%0.1f,%0.1f)" % (self.aperture.x, self.aperture.y))

        # print information about the wavelength calibration
        self.ax_header.set_xlim(0,1)
        self.ax_header.set_ylim(0,1)
        plt.setp(self.ax_header.get_yticklabels(), visible=False)
        self.ax_header.patch.set_visible(False)
        text = 'Wavelength-to-pixel guesses are {0}\n'.format(self.whichwaveid)
        text += 'Hand-picked matches are from {0}.\n'.format(self.matchesfilename.split('/')[-1])
        text += 'Pixel-to-wavelength calibration is '
        if self.justloaded:
            text += 'from ' + self.calibrationfilename.split('/')[-1]
        else:
            text += '[new!]'
        self.ax_header.text(0.025, 0.5, text,
                                    va='center', ha='left', fontsize=10)
        for i, e in enumerate(self.elements):
            self.ax_header.text(0.98,
                                1.0-(i+1.0)/(len(self.elements)+1),
                                e,
                                ha='right', va='center',
                                color=colors[e],
                                fontsize=6)

        # plot the backwards calibration
        for e in self.elements:
            ok = np.nonzero([e in self.waveids['name'][i] for i in range(len(self.waveids))])[0]

            xvals = self.waveids['pixel'][ok]
            yvals = self.waveids['wavelength'][ok]

            self.ax_header.scatter(xvals, yvals,
                            marker='o', color=colors[e], alpha=0.3)
            xfine = np.linspace(min(self.waveids['wavelength']), max(self.waveids['wavelength']))
            #self.ax_header.plot(self.waveids, xfine,
            #                    alpha=0.3)

        # plot the overlap of the lamp spectra
        scatterkw = dict(   marker='o', linewidth=0,
                            alpha=0.5, color=self.emissioncolor)
        for element in self.elements:

            self.ax_walign.plot(self.aperture.waxis,
                                self.aperture.arcs[element]['filtered'],
                                color=colors[element], alpha=0.5)
            self.ax_walign.scatter(self.pixel, self.intensity, **scatterkw)

            self.ax_walign.set_yscale('log')
            self.ax_walign.set_ylim(1, None)


        # plot tick marks for the known wavelengths
        for e in self.elements:
            for tw in self.knownwavelengths[e]:
                pix = tw['pixelguess']
                self.ax_walign.axvline(pix,
                                        ymin=0.9,
                                        color=colors[e],
                                        alpha=0.5)


        # plot the calibration
        x = np.linspace(*craftroom.oned.minmax(self.aperture.waxis),num=200)
        self.ax_wcal.plot(x, self.pixelstowavelengths(x), alpha=0.5, color='black')
        self.ax_wcal.set_ylabel('Wavelength (angstroms)')

        self.ax_wres.set_ylabel('Residuals')
        scatterkw['color'] = self.emissioncolor[self.good]
        self.ax_wcal.scatter(   self.pixel[self.good],
                                self.wavelength[self.good],
                                **scatterkw)
        self.ax_wres.scatter(   self.pixel[self.good],
                                self.residuals[self.good],
                                **scatterkw)
        scatterkw['marker'] = 'x'
        bad = self.good == False
        scatterkw['color'] = self.emissioncolor[bad]
        self.ax_wcal.scatter(self.pixel[bad], self.wavelength[bad], **scatterkw)
        self.ax_wres.scatter(self.pixel[bad], self.residuals[bad], **scatterkw)
        self.ax_wres.set_xlabel('Pixel # (by python rules)')

        self.ax_wres.set_xlim(*craftroom.oned.minmax(self.aperture.waxis))

        self.ax_walign.scatter(   self.pixel[self.handpicked],
                                self.intensity[self.handpicked],
                                marker='+', color='black')

        self.ax_wcal.scatter(   self.pixel[self.handpicked],
                                self.wavelength[self.handpicked],
                                marker='+', color='black')

        self.ax_wres.scatter(   self.pixel[self.handpicked],
                                self.residuals[self.handpicked],
                                marker='+', color='black')

        self.ax_wres.axhline(0, linestyle='--', color='gray', zorder=-100)


        rms = np.std(self.residuals[self.good])
        performance = 'using a {}-degree {}'.format(self.polynomialdegree,
                            self.pixelstowavelengths.__class__.__name__)
        performance += ' the calibration has an RMS of {0:.2f}A'.format(rms)
        performance += ' with {0:.0f} good points'.format(np.sum(self.good))
        performance += ' from {:.0f} to {:.0f}'.format(*craftroom.oned.minmax(self.pixel[self.good]))
        self.ax_wres.text(0.98, 0.05, performance,
                            fontsize=10,
                            ha='right', va='bottom',
                            transform=self.ax_wres.transAxes)

        self.ax_wres.set_ylim(*np.array([-1,1])*np.maximum(craftroom.oned.mad(self.residuals[self.good])*10, 1))
        plt.draw()
        self.speak('check out the wavelength calibration')

        if interactive == False:
            return

        options = {}
        options['q'] = dict(description='[q]uit without writing',
                            function=self.quit,
                            requiresposition=False)

        options['w'] = dict(description='[w]rite the new calibration',
                            function=self.saveandquit,
                            requiresposition=False)

        options['h'] = dict(description='match up a [h]elium line',
                            function=self.see,
                            requiresposition=True)
        options['n'] = dict(description='match up a [n]eon line',
                            function=self.see,
                            requiresposition=True)
        options['a'] = dict(description='match up a [a]rgon line',
                            function=self.see,
                            requiresposition=True)

        options['d'] = dict(description='change the polynomial [d]egree',
                            function=self.changedegree,
                            requiresposition=True)

        options['z'] = dict(description='[z]ap a previously matched line',
                            function=self.zap,
                            requiresposition=True)

        options['u'] = dict(description='[u]ndo all changes made this session',
                                    function=self.undo,
                                    requiresposition=False)

        options['r'] = dict(description='[r]estart from all defaults',
                                    function=self.restart,
                                    requiresposition=False)

        options['!'] = dict(description='raise an error[!]',
                                    function=self.freakout,
                                    requiresposition=False)

        self.speak('your options include:')
        for v in options.values():
            self.speak('   ' + v['description'])
        pressed = self.interactivewave.getKeyboard()

        try:
            # figure out which option we're on
            thing = options[pressed.key.lower()]
            # check that it's a valid position, if need be
            if thing['requiresposition']:
                assert(pressed.inaxes is not None)
            # execute the function associated with this option
            thing['function'](pressed)
        except KeyError:
            self.speak("nothing yet defined for [{}]".format(pressed.key))
            return
        except AssertionError:
            self.speak("that didn't seem to be at a valid position!")
            return

    def freakout(self, *args):
        raise RuntimeError('Breaking here, for debugginng')

    def undo(self, *args):
        self.populate()

    def restart(self, *args):
        self.populate(restart=True)

    def changedegree(self, pressed):
        '''change the degree of the polynomial'''
        self.speak('please enter a new polynomial degree (1-9)')
        new = self.interactivewave.getKeyboard()

        self.polynomialdegree = int(new.key)

    def quit(self, *args):
        '''quit without saving'''
        self.speak('finished!')
        self.notconverged = False

    def saveandquit(self, *args):
        '''save and quit'''
        self.speak('saving the wavelength calibration')
        self.save()
        self.quit()



    @property
    def handpicked(self):
        return np.array([m['mine']['handpicked'] for m in self.matches])

    def selectpeak(self, pressed):


        #print pressed.xdata, pressed.ydata
        element = shortcuts[pressed.key]
        self.speak('matching for the closest {0} line'.format(element))

        valid = np.where([e == element for e in self.pixelelement])[0]
        closest = valid[np.argmin(np.abs(self.pixel[valid] - pressed.xdata))]
        return closest


    def selectwavelength(self, pressed):
        #print(pressed.xdata, pressed.ydata)
        element = shortcuts[pressed.key.lower()]
        pixguess = [w['pixelguess'] for w in self.knownwavelengths[element]]
        closest = ((pixguess - pressed.xdata)**2).argmin()
        return self.knownwavelengths[element][closest]

    def checkvalid(self, pressed):
        if pressed.xdata is None:
            self.speak('''you typed "{0}", but the window didn't return a coordinate, please try again!'''.format(pressed.key))
            return False
        else:
            return True

    def zap(self, pressed):
        if self.checkvalid(pressed) == False:
            return
        closest = ((self.pixel - pressed.xdata)**2).argmin()
        self.matches[closest]['mine']['handpicked'] = False


    def see(self, pressed):
        if self.checkvalid(pressed) == False:
            return

        closest = self.selectpeak(pressed)
        if closest is None:
            return

        match = self.matches[closest]

        pixel = match['mine']['w']
        self.speak('you are e[X]cited about an indentifiable emission line at {0:.1}'.format(pixel))
        self.ax_walign.axvline(pixel, alpha=0.5, color=colors[match['theirs']['element']])


        self.speak('  now, please click a [H]e, [N]e, or [A]r line')

        # get a keyboard input at mouse position
        secondpressed = self.interactivewave.getKeyboard()
        if secondpressed.key.lower() in ['h', 'n', 'a']:
            # pull out the closest wavelength match
            wavematch = self.selectwavelength(secondpressed)
        else:
            self.speak("hmmm...ignoring")
            return


        self.matches[closest]['mine']['handpicked'] = True
        self.speak('updating the guess for {0}A from {1} to {2}'.format(
                        wavematch['wavelength'],
                        wavematch['pixelguess'],
                        pixel
        ))
        self.matches[closest]['theirs']['pixelguess'] = pixel
        self.matches[closest]['mine']['wavelength'] = wavematch['wavelength']


    def exclude(self, pressed):
        self.speak('excluding a point')
        if pressed.inaxes == self.ax_wres:
            closest = self.select(pressed)
            self.matches[closest]['mine']['handpicked'] = False
