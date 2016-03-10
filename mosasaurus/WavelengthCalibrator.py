'''Wavelength Calibrator defines the wavelength calibration, and '''

from imports import *
colors = dict(He='lightsalmon', Ne='red', Ar='deepskyblue')
shortcuts = {'h':'He', 'n':'Ne', 'a':'Ar'}

class WavelengthCalibrator(Talker):
    def __init__(self, aperture):
        Talker.__init__(self)
        self.aperture = aperture
        self.filename = self.aperture.directory + 'waveCal_{0}.npy'.format(self.aperture.name)

        self.waveids =  astropy.io.ascii.read(self.aperture.obs.wavelength2pixelsFile)

        self.elements = ['He', 'Ne','Ar']

        self.wavelengthorder = 3
        self.matchdistance = 20
        self.load()

    def save(self):
        np.save(filename, self.waveCalCoef)
        self.speak("saved wavelength calibration to {0}".format(self.filename))

    def load(self):
        try:
            self.coef = np.load(self.filename)
            self.speak("loaded wavelength calibration from {0}".format(self.filename))
            self.justloaded = True

            # need the inverse relation for making the plots make sense
            inversecoef = np.polyfit(self.pixelstowavelengths(self.aperture.waxis), self.aperture.waxis, self.wavelengthorder)
            self.wavelengthstopixels = np.poly1d(inversecoef)

        except IOError:
            self.justloaded = False
        self.create()


    def findRoughShift(self, blob=2.0):
        '''using a list of pixels matched to wavelengths, find the rough
                offset of the this arc relative to the standard slits/setting'''

        self.speak("cross correlating arcs with known wavelengths")

        # create a plot showing how well the lines match
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
            x = np.arange(-self.aperture.obs.ysize,self.aperture.obs.ysize)

            myPeaks, theirPeaks = np.zeros(len(x)), np.zeros(len(x))
            # create fake spectrum of their peaks
            for i in range(len(self.waveids)):
                if element in self.waveids['name'][i]:
                    center = self.waveids['pixels'][i]/self.aperture.obs.binning
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

        # find the peak of the combined correlation function
        self.peakoffset = np.where(self.corre['combined'] == self.corre['combined'].max())[0][0] - len(x)
        # (old?) to convert: len(x) - xPeak = x + peakoffset

        # plot the shifted wavelength ids, and combined corfuncs
        for element in self.elements:
            for i in range(len(self.waveids)):
                if element in self.waveids['name'][i]:
                      center = self.waveids['pixels'][i]/self.aperture.obs.binning
                      center += self.peakoffset
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
            self.ax_wavecor[self.elements[0]].set_title('cross correlation peaks at \n{0} pixels ({1}x{1} binned pixels)'.format(self.peakoffset, self.aperture.obs.binning), fontsize=fontsize)
            self.ax_waverough[self.elements[0]].set_title(
            'Coarse Wavelength Alignment\nfor ({0:0.1f},{1:0.1f})'.format(
                 self.aperture.x, self.aperture.y),fontsize=fontsize)

            # save the figure
            figure_waverough.savefig(
            self.aperture.directory + 'roughWavelengthAlignment_{0}.pdf'.format(
            self.aperture.name))

    def findPeaks(self):
        '''identify peaks in the extracted arc spectrum'''

        # extract a spectrum from the master image for each lamp
        self.aperture.arcs = {}
        for element in self.elements:
            self.aperture.arcs[element] = self.aperture.extract(
                                image=self.aperture.images[element],
                                arc=True)

        # load the complete list of wavelengths (with no pixels)
        self.waveall=astropy.io.ascii.read(self.aperture.obs.wavelengthsFile)

        # find the peaks in my spectra
        self.peaks = {}

        for count, element in enumerate(self.elements):

            # the pixel spectrum self.aperture.extracted from this arc lamp
            flux = self.aperture.arcs[element]['raw_counts']

             # identify my peaks
            xPeak, yPeak = zachopy.oned.peaks(self.aperture.waxis, flux)

            # for some reason, need to trim peaks outside range
            pad = 25
            toright = xPeak > (np.min(self.aperture.waxis) + pad)
            toleft = xPeak < (np.max(self.aperture.waxis) - pad)
            ok = toleft*toright
            xPeak, yPeak = xPeak[ok], yPeak[ok]

            n = len(xPeak)

            # store those peaks
            self.peaks[element] = []
            for i in range(n):
                peak =  {
                        'w':xPeak[i],
                        'intensity':yPeak[i],
                        'handpicked':False,
                        'outlier':False
                        }
                self.peaks[element].append(peak)

    def findKnownWavelengths(self, extrapolateorder=3):
        # create a temporary calibration to match reference wavelengths to reference pixels (so we can extrapolate to additional wavelengths not recorded in the dispersion solution file)

        coef = np.polyfit(
            (self.waveids['wave2']),
            self.waveids['pixels']/self.aperture.obs.binning+self.peakoffset,
            extrapolateorder)
        self.wavelengthstopixels = np.poly1d(coef)


        self.knownwavelengths = {}

        # treat the arc lamps separately
        for count, element in enumerate(self.elements):
            # pull out the wavelengths from the complete file
            self.knownwavelengths[element] = []
            for i in range(len(self.waveall)):
                if element in self.waveall['name'][i]:
                    wave = self.waveall['wavelength'][i]
                    known = {
                            'element':element,
                            'wavelength':wave,
                            'pixelguess':self.wavelengthstopixels(wave)
                            }
                    self.knownwavelengths[element].append(known)

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
                closest = np.nonzero(np.abs(distance) == np.min(np.abs(distance)))[0]

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
                    print "{0} of my peaks match to their {1}".format(len(relevant), theirs)

                    distances = np.abs([m['theirs']['pixelguess'] - m['mine']['w'] for m in relevant])


                    best = distances.argmin()

                    print "the closests one is {0}".format(relevant[best]['mine']['w'])
                    self.matches.append(relevant[best])
            '''



            # add this to the list
            #self.matches.extend(thiselement)

    def create(self, remake=False):
        '''Populate the wavelength calibration for this aperture.'''

        self.speak("populating wavelength calibration")

        # pull out the peaks from extracted arc spectra
        self.findPeaks()

        # perform a first rough alignment, to make pixels
        self.findRoughShift()

        self.findKnownWavelengths()

        # loop through
        self.notconverged = True
        while(self.notconverged):

            # set an initial guess matching wavelengths to my pixles
            self.guessMatches()

            # do a fit
            if self.justloaded:
                self.justloaded = False
            else:
                # do an initial fit
                self.coef = np.polyfit(     self.pixel,
                                            self.wavelength,
                                            self.wavelengthorder,
                                            w=self.weights)
            # identify outliers
            limit = 1.48*zachopy.oned.mad(self.residuals[self.good])*4
            limit = np.maximum(limit, 1.0)
            # outliers get reset each time (so don't lose at edges)
            outlier = np.abs(self.residuals) > limit
            # keep track of which are outlier
            for i, m in enumerate(self.matches):
                self.matches[i]['mine']['outlier'] = outlier[i]
            self.speak('points beyond {0} are outliers ({1})'.format(limit, i))
            self.coef = np.polyfit(     self.pixel,
                                        self.wavelength,
                                        self.wavelengthorder,
                                        w=self.weights)

            self.plotWavelengthFit()
            self.updatew2p()

    @property
    def weights(self):
        return self.good+self.handpicked*10

    def updatew2p(self):
        coef = np.polyfit(self.wavelength, self.pixel,
                            self.wavelengthorder, w=self.weights)
        self.wavelengthstopixels = np.poly1d(coef)
        #self.matchdistance = 2*np.std((self.pixel - self.wavelengthstopixels(self.wavelength))[self.good])
        #self.matchdistance = np.maximum(self.matchdistance, 10)
        self.speak('match distance is {0}'.format(self.matchdistance))

    @property
    def residuals(self):
        return  self.wavelength - self.pixelstowavelengths(self.pixel)

    @property
    def good(self):
        isntoutlier = np.array([m['mine']['outlier'] == False for m in self.matches])
        ishandpicked = np.array([m['mine']['handpicked'] for m in self.matches])
        return (isntoutlier + ishandpicked) > 0

    @property
    def pixelstowavelengths(self):
        return np.poly1d(self.coef)

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

    def plotWavelengthFit(self):
        # plot to make sure the wavelength calibration makes sense

        figure_wavelengthcal = plt.figure('wavelength calibration',
                                        figsize=(10,4), dpi=100)
        self.interactivewave = zachopy.iplot.iplot(4,1,
                height_ratios=[0.1, 0.4, 0.2, .2], hspace=0)

        self.ax_w2p = self.interactivewave.subplot(0)

        # do the lamp spectra overlap?
        self.ax_walign = self.interactivewave.subplot(1, sharex=self.ax_w2p)

        # what is the actual wavelength calibration
        self.ax_wcal = self.interactivewave.subplot(2, sharex=self.ax_w2p)

        # what are the residuals from the fit
        self.ax_wres = self.interactivewave.subplot(3, sharex=self.ax_w2p)

        for ax in [self.ax_w2p, self.ax_walign, self.ax_wcal]:
            plt.setp(ax.get_xticklabels(), visible=False)

        self.ax_w2p.set_title("Wavelength Calib. for Aperture"
            " (%0.1f,%0.1f)" % (self.aperture.x, self.aperture.y))

        kw = dict(marker='o')

        # plot the backwards calibration
        for e in self.elements:
            ok = np.nonzero([e in self.waveids['name'][i] for i in range(len(self.waveids))])[0]

            xvals = self.waveids['pixels'][ok]/self.aperture.obs.binning + self.peakoffset
            yvals = self.waveids['wave2'][ok]

            self.ax_w2p.scatter(xvals, yvals,
                            marker='o', color=colors[e], alpha=0.3)
            xfine = np.linspace(min(self.waveids['wave2']), max(self.waveids['wave2']))
            self.ax_w2p.plot(self.wavelengthstopixels(xfine), xfine,
                                alpha=0.3)

        # plot the overlap of the lamp spectra
        for element in self.elements:

            self.ax_walign.plot(self.aperture.waxis,
                                self.aperture.arcs[element]['raw_counts'],
                                color=colors[element], alpha=0.5)
            self.ax_walign.scatter(self.pixel, self.intensity,
                                    marker='o', linewidth=0,
                                    alpha=0.5, color=self.emissioncolor)

            self.ax_walign.set_yscale('log')


        # plot tick marks for the known wavelengths
        for e in self.elements:
            for tw in self.knownwavelengths[e]:
                pix = self.wavelengthstopixels(tw['wavelength'])
                self.ax_walign.axvline(pix,
                                        ymin=0.9,
                                        color=colors[e],
                                        alpha=0.5)

        # plot the new calibration
        self.ax_wcal.scatter(self.pixel, self.wavelength, color=self.emissioncolor, **kw)
        self.ax_wcal.plot(self.pixel, self.pixelstowavelengths(self.pixel), alpha=0.5, color='black')
        self.ax_wcal.set_ylabel('Wavelength (angstroms)')

        self.ax_wres.set_ylabel('Residuals')
        self.ax_wres.scatter(self.pixel[self.good], self.residuals[self.good], color=self.emissioncolor[self.good], **kw)
        kw['marker'] = 'x'
        bad = self.good == False
        self.ax_wres.scatter(self.pixel[bad], self.residuals[bad], color=self.emissioncolor[bad], **kw)
        self.ax_wres.set_xlabel('Pixel # (by python rules)')

        self.ax_wres.set_xlim(min(self.aperture.waxis), max(self.aperture.waxis))

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
        self.ax_wres.set_ylim(*np.array([-1,1])*np.maximum(zachopy.oned.mad(self.residuals[self.good])*10, 1))
        plt.draw()
        self.speak('check out the wavelength calibration')
        pressed = self.interactivewave.getKeyboard()

        if pressed.key == 'q':
            # quit, deciding this was great
            self.notconverged = False
            figure_wavelengthcal.savefig(self.aperture.directory + 'wavelengthCalibration_{0}.pdf'.format(self.aperture.name))
        elif pressed.key in ['h', 'n', 'a']:
            # see that a line could be made to line up
            self.see(pressed)
        elif pressed.key == 'o':
            self.wavelengthorder = int(self.input('please enter a new order!'))
        elif pressed.key == 'z':
            # zap a previously pinned point
            self.zap(pressed)
        elif pressed.key == 'r':
            # zap a previously pinned point
            self.create()
        elif pressed.key == '!':
            # quit, and raise an error
            assert(False)
        else:
            return

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
        #print pressed.xdata, pressed.ydata
        element = shortcuts[pressed.key.lower()]
        pixguess = [w['pixelguess'] for w in self.knownwavelengths[element]]
        closest = ((pixguess - pressed.xdata)**2).argmin()
        return self.knownwavelengths[element][closest]


    def zap(self, pressed):
        closest = ((self.pixel - pressed.xdata)**2).argmin()
        self.matches[closest]['mine']['handpicked'] = False


    def see(self, pressed):
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
