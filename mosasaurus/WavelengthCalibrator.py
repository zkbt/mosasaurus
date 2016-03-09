'''Wavelength Calibrator defines the wavelength calibration, and '''

from imports import *
colors = dict(He='lightsalmon', Ne='red', Ar='deepskyblue')


class WavelengthCalibrator(Talker):
    def __init__(self, aperture):
        Talker.__init__(self)
        self.aperture = aperture
        self.waveids =  astropy.io.ascii.read(self.aperture.obs.wavelength2pixelsFile)

        self.elements = ['He', 'Ne','Ar']
        self.wavelengthorder = 3
        self.matchdistance = 20
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
        ax_waverough, ax_wavecor = {}, {}
        sharer, sharec = None, None
        for i, e in enumerate(self.elements):
            ax_waverough[e] = plt.subplot(gs[i,0],sharex=sharer)
            ax_wavecor[e] = plt.subplot(gs[i,1], sharex=sharec)
            sharer, sharec = ax_waverough[e], ax_wavecor[e]

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
            ax_waverough[element].plot(x, myPeaks/myPeaks.max(),
                                            label='extracted', alpha=0.5,
                                            color=colors[element])
            ax_waverough[element].set_ylim(0, 2)


            # plot the correlation functions
            normcor = self.corre[element]
            assert(np.isfinite(self.corre[element]).any())

            normcor /= np.nanmax(self.corre[element])
            ax_wavecor[element].plot(normcor, label=element, alpha=0.5,
                                            color=colors[element])
            # tidy up the plots
            for a in [ax_wavecor, ax_waverough]:
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
                          ax_waverough[element].axvline(center, alpha=0.25, color='black')
                # plot the combined correlation function
                normedcombined = self.corre['combined']/np.max(self.corre['combined'])
                ax_wavecor[element].plot(normedcombined,
                                    label='combined', alpha=0.25, color='black')
                # tidy up the plots
                ax_waverough[element].set_ylabel(element)
                ax = ax_waverough[self.elements[-1]]
            plt.setp(ax.get_xticklabels(), visible=True)
            ax.set_xlabel('Pixel Position')
            fontsize = 8
            ax_wavecor[self.elements[0]].set_title('cross correlation peaks at \n{0} pixels ({1}x{1} binned pixels)'.format(self.peakoffset, self.aperture.obs.binning), fontsize=fontsize)
            ax_waverough[self.elements[0]].set_title(
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

            # loop over my peaks
            for i in range(len(theirPeaksOnMyPixels)):

                # find my closest peak to theirs
                distance = myPeaks - theirPeaksOnMyPixels[i]
                closest = np.nonzero(np.abs(distance) == np.min(np.abs(distance)))[0]

                if distance[closest] < self.matchdistance:
                    # call this a match (no matter how far)
                    match = {
                            'mine':self.peaks[element][closest],
                            'theirs':self.knownwavelengths[element][i],
                            'distance':distance[closest]
                            }

                    # add this to the list
                    self.matches.append(match)

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
            self.coef = np.polyfit(     self.pixel,
                                        self.wavelength,
                                        self.wavelengthorder,
                                        w=self.weights)
            print self.good
            # identify outliers
            limit = 1.48*zachopy.oned.mad(self.residuals[self.good])*2
            limit = np.maximum(limit, 1.0)
            outlier = np.abs(self.residuals) > limit
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
    def wavelength(self):
        return np.array([m['theirs']['wavelength'] for m in self.matches])

    @property
    def emissioncolor(self):
        return np.array([colors[m['theirs']['element']] for m in self.matches])

    def plotWavelengthFit(self):
        # plot to make sure the wavelength calibration makes sense

        figure_wavelengthcal = plt.figure('wavelength calibration',
                                        figsize=(6,4), dpi=100)
        interactivewave = zachopy.iplot.iplot(4,1,
                height_ratios=[.2, .2, .2, .4], hspace=0)

        # does the wavelength2pixel code work
        ax_w2p = interactivewave.subplot(0)

        # do the lamp spectra overlap?
        ax_walign = interactivewave.subplot(1, sharex=ax_w2p)

        # what is the actual wavelength calibration
        ax_wcal = interactivewave.subplot(2, sharex=ax_w2p)

        # what are the residuals from the fit
        self.ax_wres = interactivewave.subplot(3, sharex=ax_w2p)

        for ax in [ax_w2p, ax_walign, ax_wcal]:
            plt.setp(ax.get_xticklabels(), visible=False)

        ax_w2p.set_title("Wavelength Calib. for Aperture"
            " (%0.1f,%0.1f)" % (self.aperture.x, self.aperture.y))

        kw = dict(marker='o')

        # plot the backwards calibration
        for e in self.elements:
            ok = np.nonzero([e in self.waveids['name'][i] for i in range(len(self.waveids))])[0]

            xvals = self.waveids['pixels'][ok]/self.aperture.obs.binning + self.peakoffset
            yvals = self.waveids['wave2'][ok]

            ax_w2p.scatter(xvals, yvals, marker='o', color=colors[e])
            xfine = np.linspace(min(self.waveids['wave2']), max(self.waveids['wave2']))
            ax_w2p.plot(self.wavelengthstopixels(xfine), xfine)

        # plot the overlap of the lamp spectra
        for element in self.elements:

            ax_walign.plot(self.aperture.waxis, self.aperture.arcs[element]['raw_counts'],
                                color=colors[element], alpha=0.5)
            ax_walign.set_yscale('log')


        for i,tw in enumerate(self.waveall['wavelength']):
            name = self.waveall['name'][i][0:2]
            if name in self.elements:
                ax_walign.axvline(self.wavelengthstopixels(tw), ymin=0.9,
                                        color=colors[name],
                                        alpha=0.5)

        # plot the new calibration
        ax_wcal.scatter(self.pixel, self.wavelength, color=self.emissioncolor, **kw)
        ax_wcal.plot(self.pixel, self.pixelstowavelengths(self.pixel), alpha=0.5, color='black')
        ax_wcal.set_ylabel('Wavelength (angstroms)')

        self.ax_wres.set_ylabel('Residuals')
        self.ax_wres.scatter(self.pixel[self.good], self.residuals[self.good], color=self.emissioncolor[self.good], **kw)
        kw['marker'] = 'x'
        bad = self.good == False
        self.ax_wres.scatter(self.pixel[bad], self.residuals[bad], color=self.emissioncolor[bad], **kw)
        self.ax_wres.set_xlabel('Pixel # (by python rules)')
        self.ax_wres.set_xlim(min(self.aperture.waxis), max(self.aperture.waxis))

        self.ax_wres.scatter(   self.pixel[self.handpicked],
                                self.residuals[self.handpicked],
                                marker='+', color='black')

        self.ax_wres.axhline(0, linestyle='--', color='gray', zorder=-100)
        self.ax_wres.set_ylim(*np.array([-1,1])*zachopy.oned.mad(self.residuals)*10)
        plt.draw()
        self.speak('check out the wavelength calibration')
        pressed = interactivewave.getKeyboard()

        if pressed.key == 'q':
            self.notconverged = False
            figure_wavelengthcal.savefig(self.aperture.directory + 'wavelengthCalibration_{0}.pdf'.format(self.aperture.name))
        elif pressed.key == 'c':
            self.choose(pressed)
        elif pressed.key == 'x':
            self.exclude(pressed)
        elif pressed.key == '!':
            assert(False)
        else:
            return

    @property
    def handpicked(self):
        return np.array([m['mine']['handpicked'] for m in self.matches])

    def select(self, pressed):
        #print pressed.xdata, pressed.ydata
        return (((self.pixel - pressed.xdata)**2 + (self.residuals - pressed.ydata))**2).argmin()

    def choose(self, pressed):
        self.speak('choosing a point')
        if pressed.inaxes == self.ax_wres:
            closest = self.select(pressed)
            self.matches[closest]['mine']['handpicked'] = True

    def exclude(self, pressed):
        self.speak('excluding a point')
        if pressed.inaxes == self.ax_wres:
            closest = self.select(pressed)
            self.matches[closest]['mine']['handpicked'] = False
