'''Trace defines the extraction and sky subtraction regions for an Aperture'''

from .imports import *
from craftroom.cmaps import one2another

class Trace(Talker):
    '''trace object defines the extraction trace and sky region'''
    def __init__(self, aperture):
        Talker.__init__(self)
        self.aperture = aperture
        self.obs = self.aperture.obs
        self.instrument = self.obs.instrument

        # keep track of the spatial and wavelength axis (in pixels)
        self.saxis = self.aperture.saxis
        self.waxis = self.aperture.waxis

        # what degree polynomial should be used?
        self.order = self.instrument.extractiondefaults['traceOrder']

        # define a grid of extraction widths, (can be changed later, if desired)
        self.setSizes(default=True)
        self.traceWidth = np.min(self.extractionwidths)

        # keep track of where the crosshair is pointing
        self.crosshair = dict(s=0.0, w=0.0)

        try:
            self.speak('trying to load trace parameters from {}'.format(self.filename))
            self.load()
        except IOError:
            self.speak('interactively creating a new trace')

            # set up the initial sky offsets
            inner = np.min(self.extractionwidths) + self.instrument.extractiondefaults['skyGap']
            outer = inner + self.instrument.extractiondefaults['skyWidth']
            self.skyoffsets = [ dict(top=outer, bottom=inner, whattodo=1),
                                dict(top=-inner, bottom=-outer, whattodo=1)]


            #maskinner = outer + self.instrument.extractiondefaults['skyGap']
            #maskouter = maskinner + self.instrument.extractiondefaults['skyWidth']
            #self.maskfluxoffsets = [ dict(top=maskouter, bottom=maskinner, whattodo=1),
            #                         dict(top=-maskinner, bottom=-maskouter, whattodo=1)]

            # an array of custom-selected trace points (for forcing a particular fit)
            self.tracepoints = {'w':[], 's':[]}

            # and the initial guess for the trace
            self.traceguess = np.poly1d(0.0)

            # set up the intial bits of the trace object
            self.setup()

            # run the interactive loop
            self.run()

        # save the defined properties of the trace
        # self.save() # (don't need, because user will be quitting with "w"?)

    def writeandquit(self, *args):
        '''write out all the parameters to a file, and stop the interaction loop'''
        # save the parameters
        self.save()
        self.quit()

    def quit(self, *args):
        self.notconverged = False
        plt.close(self.figure)

    @property
    def filename(self):
        return os.path.join(self.aperture.directory,  'trace_{0}.npy'.format(self.aperture.name))

    def save(self):
        '''save all properties of this trace to a file'''

        # save the parameters of the trace
        np.save(self.filename, (self.tracefitcoef, self.tracefitwidth))
        self.speak("saved trace parameters to {0}".format(self.filename))

        # save the extraction and sky masks
        filename = os.path.join(self.aperture.directory, 'extractionmasks_{0}.npy'.format(self.aperture.name))
        np.save(filename, (self.skyoffsets, self.extractionwidths))
        #np.save(filename, (self.skyoffsets, self.maskfluxoffsets, self.extractionwidths))
        self.speak("saved extraction mask parameters to {0}".format(filename))

        # save a PDF of the trace definition
        filename = os.path.join(self.aperture.directory, 'tracedefinition_{0}.pdf'.format(self.aperture.name))
        self.figure.savefig(filename)

        # save skyMask and extractionMask for each width
        filename = os.path.join(self.aperture.directory, 'masks_{0}.npy'.format(self.aperture.name))
        masks = {}
        for width in self.extractionwidths:
            masks[width] = {}
            masks[width]['extractionMask'] =  self.extractionmask(width)
            masks[width]['skyMask'] = self.skymask(width)
        np.save(filename, masks)
        self.speak("saved set of extraction and sky masks for each width")            

    def load(self):

        # load the parameters of the trace
        (self.tracefitcoef, self.tracefitwidth) = np.load(self.filename, allow_pickle=True)
        self.speak("loaded trace parameters to {0}".format(self.filename))
        self.tracefit = np.poly1d(self.tracefitcoef)

        # save the extraction and sky masks
        filename = os.path.join(self.aperture.directory,  'extractionmasks_{0}.npy'.format(self.aperture.name))
        #(self.skyoffsets, self.maskfluxoffsets, self.extractionwidths) = np.load(filename)
        (self.skyoffsets, self.extractionwidths) = np.load(filename, allow_pickle=True)
        self.speak("saved extraction mask parameters to {0}".format(filename))
        self.numberofapertures = len(self.extractionwidths)
        self.narrowest, self.widest = craftroom.oned.minmax(self.extractionwidths)

    def run(self):
        '''interactively fit for the position of the trace of the spectrum,
                its width, and some sky subtraction regions

                when running, the following commands can be used:

                    press "t" at any point on the spectral trace

                        Add point to the object's trace.
                        It doesn't need to be exact, this
                        is simply to constrain an initial
                        guess for the trace. Centroids and
                        widths will be calculated to determine
                        the optimal trace coefficients.

                    press "e" at a point on the sky,
                    and again at a different point on the sky

                        "Extend" the sky subtraction region.
                        All this cares about is the offset
                        in the cross-dispersion direction
                        away from the spectral trace.

                    press "r" at a point on the sky,
                    and again at a different point on the sky

                        "Remove" from the sky subtraction region.
                        All this cares about is the offset
                        in the cross-dispersion direction
                        away from the spectral trace.

                    press "c" at a point on the image

                        Move the "crosshair" to an (x,y)
                        position, and plot the corresponding
                        slices on the projected axes.

                    press "f" anywhere in the image region

                        "Fit" for the spectral trace, and
                        sky subtraction regions, using the
                        current inputs.

                    (something to set the width[s])

                    (something to change the plotting scale)
                '''
        self.notconverged = True
        while self.notconverged:
            self.speak('Please refine the extraction aperture.')

            options = {}
            options['w'] = dict(description="[w]rite and quit",
                                function=self.writeandquit,
                                requiresposition=False)
            options['q'] = dict(description="[q]uit without writing",
                                function=self.quit,
                                requiresposition=False)
            options['c'] = dict(description="move the [c]rosshair, and plot slices along it",
                                function=self.moveCrosshair,
                                requiresposition=True)
            options['t'] = dict(description="add a guess for a [t]race point",
                                function=self.addTracePoint,
                                requiresposition=True)
            options['e'] = dict(description="[e]xtend a sky region (twice for start and stop)",
                                function=self.modifySky,
                                requiresposition=True)
            options['r'] = dict(description="[r]emove a sky region (twice for start and stop)",
                                function=self.modifySky,
                                requiresposition=True)
            options['f'] = dict(description="[f]it the trace, using the star's centroids and sky areas",
                                function=self.fitTrace,
                                requiresposition=False)
            options['s'] = dict(description="[s]et the [s]ize of the [s]tar's smallest extraction region",
                                function=self.setSizes,
                                requiresposition=False)

            # hzdl:
            # trying to add designation for "out of slit" or "mask" flux.
            # want to designate this similar to sky background
            #options['m'] = dict(description='extend a [m]ask region (twice for start and stop)',
            #                    function=self.modifyMaskFlux,
            #                    requiresposition=True)
            #options['n'] = dict(description="[n]o, not that mask region, like 'r' (twice for start and stop)",
            #                    function=self.modifyMaskFlux,


            # print the options
            self.speak('your options include:')
            for v in options.values():
                self.speak('   ' + v['description'])

            # get the keyboard input
            pressed = self.iplot.getKeyboard()

            # process the keyboard input
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
            except AssertionError:
                self.speak("that didn't seem to be at a valid position!")


            plt.draw()

    @property
    def traceCenter(self):
        '''return a function that gives the s-value of the trace center, given an input w'''
        try:
            return self.tracefit
        except AttributeError:
            return self.traceguess

    def setSizes(self, pressed=None, default=False):
        '''prompt the user to select range of aperture sizes for extraction'''
        if default:
            self.narrowest = self.instrument.extractiondefaults['narrowest']
            self.widest = self.instrument.extractiondefaults['widest']
            self.numberofapertures = self.instrument.extractiondefaults['numberofapertures']
        else:
            self.speak('Please redefine the aperture sizes.')
            self.speak(' (The old aperture sizes were {}.)'.format(self.extractionwidths))
            self.speak(' ["d" for any of these for defaults]')

            try:
                self.narrowest = np.float(self.input("What is the narrowest aperture?"))
            except ValueError:
                self.narrowest = self.instrument.extractiondefaults['narrowest']

            try:
                self.widest = np.float(self.input("What is the widest aperture?"))
            except ValueError:
                self.widest = self.instrument.extractiondefaults['widest']

            try:
                self.numberofapertures = np.int(self.input('How many apertures do you want?'))
            except ValueError:
                self.numberofapertures = self.instrument.extractiondefaults['numberofapertures']

        self.extractionwidths = np.linspace(self.narrowest, self.widest, self.numberofapertures)
        self.speak(' The current aperture sizes are {}.'.format(self.extractionwidths))

        # update the plotting, to reflect the new extraction apertures
        try:
            self.updateMasks()
        except AttributeError:
            pass

    def modifySky(self, pressed, whattodo=None):
        '''from a KeyEvent, extend or remove regions from the sky'''
        try:
            # if we've clicked once before, add to those positions
            self.skykeys.append(pressed)
        except AttributeError:
            # if we haven't clicked yet, create a new catalog
            self.skykeys = [pressed]

        # pull out the values
        w = np.array([sk.xdata for sk in self.skykeys])
        s = np.array([sk.ydata for sk in self.skykeys])
        self.plotted['dragger'].set_data(w, s)

        if len(self.skykeys) == 2:
            a, b = [sk.key for sk in self.skykeys]
            if a != b:
                self.speak("Uh-oh, you pressed {0} and then {1}, "
                    "but that doesn't make sense. "
                    "Please try again!".format(a,b))
                return

            # are we extending or removing?
            whattodo = {"e":1, "r":0}[a]

            # figure out offset from the trace guess
            offsets = s - self.traceCenter(w)

            d = {   'top':offsets.max(),
                    'bottom':offsets.min(),
                    'whattodo':whattodo}

            try:
                self.skyoffsets.append(d)
            except AttributeError:
                self.skyoffsets = [d]

            self.updateMasks()

            # remove the dictionary, to start over again for the next
            del self.skykeys

        elif len(self.skykeys) == 1:
            self.speak('Please do another "{}" again to edit the sky.'.format(pressed.key))
        else:
            self.speak('How did you get here?')

    #def modifyMaskFlux(self, pressed, whattodo=None):
    #    '''from a KeyEvent, extend or remove regions from the mask to be used to estimate extra out-of-slit flux '''
    #    try:
    #        # if we've clicked once before, add to those positions
    #        self.maskfluxkeys.append(pressed)
    #    except AttributeError:
    #        # if we haven't clicked yet, create a new catalog
    #        self.maskfluxkeys = [pressed]

        # pull out the values
    #    w = np.array([mf.xdata for mf in self.maskfluxkeys])
    #    s = np.array([mf.ydata for mf in self.maskfluxkeys])
    #    self.plotted['dragger'].set_data(w, s)

    #    if len(self.maskfluxkeys) == 2:
    #        a, b = [mf.key for mf in self.maskfluxkeys]
    #        if a != b:
    #            self.speak("Uh-oh, you pressed {0} and then {1}, "
    #                "but that doesn't make sense. "
    #                "Please try again!".format(a,b))
    #            return

            # are we extending or removing?
    #        whattodo = {"m":1, "n":0}[a]

            # figure out offset from the trace guess
    #        offsets = s - self.traceCenter(w)

    #        d = {   'top':offsets.max(),
    #                'bottom':offsets.min(),
    #                'whattodo':whattodo}

    #        try:
    #            self.maskfluxoffsets.append(d)
    #        except AttributeError:
    #            self.maskfluxoffsets = [d]

    #        self.updateMasks()

            # remove the dictionary, to start over again for the next
    #        del self.maskfluxkeys

    #    elif len(self.maskfluxkeys) == 1:
    #        self.speak('Please do another "{}" again to edit the mask flux.'.format(pressed.key))
    #    else:
    #        self.speak('How did you get here?')

    def updateMasks(self):
        '''update the sky mask (because either the sky offsets have changed),
            or the trace center has changed), and update the plots'''


        # update the extraction edges
        for line in self.plotted['extractionedges']:
            line.remove()
        offsets = (np.array([-1,1])[np.newaxis, :]*self.extractionwidths[:,np.newaxis]).flatten()
        x = self.waxis
        y = self.traceCenter(x)[:,np.newaxis] + offsets[np.newaxis, :]
        self.plotted['extractionedges'] = self.ax['2d'].plot(x, y,
                                            linewidth=1,
                                            alpha=0.5,
                                            color='darkgreen',
                                            zorder = 1000)

        self.plotted['extractionmask'].set_data(self.extractionmaskimage)
        self.plotted['skymask'].set_data(self.skymaskimage)
        #self.plotted['maskfluxmask'].set_data(self.maskfluxmaskimage)

        for t in self.plotted['widthlabels']:
            t.remove()


        self.labelWidths()
        plt.draw()

    def addTracePoint(self, pressed):
        '''from a KeyEvent, add a point to the list of trace guesses.
        this will also update the trace to fit the guesses, and not the centroids.
        run "f" again to snap back to the centroiding fit'''

        x,y = pressed.xdata, pressed.ydata
        self.speak('adding trace guess point at {0}'.format((x,y)))

        # add the trace points to the list of guesses
        self.tracepoints['w'].append(x)
        self.tracepoints['s'].append(y)

        # update the plot
        self.plotted['tracepoints'].set_data(self.tracepoints['w'],
                                             self.tracepoints['s'])

        # update the guess
        order = np.minimum(len(self.tracepoints['w'])-1, self.order)
        self.traceguesscoef = np.polyfit(self.tracepoints['w'], self.tracepoints['s'], order)
        #self.speak('the trace guess coefficients are {}'.format(self.traceguesscoef))
        self.traceguess = np.poly1d(self.traceguesscoef)

        xfine = np.linspace(self.waxis.min(), self.waxis.max(), 10)
        self.plotted['traceguess'].set_data(xfine, self.traceguess(xfine))

        # set the trace fit to be the guess!
        self.tracefit = self.traceguess
        self.moveCrosshair(pressed)

        # the sky moves, if the trace does
        self.updateMasks()

    def fitTrace(self, *args):
        '''using the guess of the trace, and the sky subtraction regions,
            fit for trace using actual centroids'''

        # estimate a rough 1D spectrum
        flattened = self.aperture.images['science']/self.aperture.images['RoughFlat']
        roughSky1d=np.average(flattened,
                                axis=self.aperture.sindex,
                                weights=self.skymask(np.median(self.extractionwidths)))
        reshaped = roughSky1d.reshape((self.waxis.shape[0],1))
        self.aperture.images['Sky'] = reshaped*np.ones_like(flattened)
        self.aperture.images['Subtracted'] = flattened - self.aperture.images['Sky']
        fine = np.isfinite(self.aperture.images['Subtracted'])


        considerstar = self.extractionmask(np.max(self.extractionwidths))

        weights = np.maximum(fine*considerstar*self.aperture.images['Subtracted'], 0)
        weights += 1.0/np.sum(weights)


        # use the roughly subtracted image to fit centroids
        flux = np.average(self.aperture.images['Subtracted'],
                            axis=self.aperture.sindex,
                            weights=weights)

        if (np.sum(weights, self.aperture.sindex) <= 0).any():
            self.speak('weights were wonky')
            assert(False)

        try:
            fluxWeightedCentroids = np.ma.average(self.aperture.s,
                        axis=self.aperture.sindex,
                        weights=fine*considerstar*self.aperture.images['Subtracted'])
        except ZeroDivisionError:
            self.speak('UH-OH, got zero-division error')
            return
        if np.isfinite(fluxWeightedCentroids).all() == False:
            self.speak('centroids were wacky')
            assert(False)
        reshapedFWC = fluxWeightedCentroids[:,np.newaxis]
        weights = fine*considerstar*self.aperture.images['Subtracted']
        weights = np.maximum(weights, 0)
        weights += 1.0/np.sum(weights)
        assert(weights.sum() > 0)
        fluxWeightedWidths = np.sqrt(
                    np.average((self.aperture.s-reshapedFWC)**2, axis=self.aperture.sindex, weights=weights))
        assert(np.isfinite(fluxWeightedWidths).all())

        fluxWeightedWidths[np.isfinite(fluxWeightedWidths)==False] = np.inf

        self.plotted['fittedwidthsup'].set_data(self.waxis,
                    fluxWeightedCentroids + fluxWeightedWidths)
        self.plotted['fittedwidthsdown'].set_data(self.waxis,
                    fluxWeightedCentroids - fluxWeightedWidths)

        # fit a polynomial to the ridge of the spectrum
        self.tracefitcoef = np.polyfit(self.waxis, fluxWeightedCentroids,
                            self.order, w=flux/fluxWeightedWidths**2)


        self.tracefit = np.poly1d(self.tracefitcoef)
        self.tracefitwidth = np.median(fluxWeightedWidths)
        self.traceWidth = self.tracefitwidth*3
        self.updateMasks()

    def setup(self, percentiles=[1,99]):
        '''make the plotting window and interactive tools'''

        plt.ion()
        self.figure = plt.figure('tracing the spectrum',
                                    figsize=(8,4), dpi=100)

        # create an interactive plot
        self.iplot = craftroom.displays.iplot.iplot(2,2,
                                        hspace=0, wspace=0,
                                        left=0.05, right=0.95,
                                        bottom=0.05, top=0.95,
                                        width_ratios=[1.0, 0.1],
                                        height_ratios=[0.1, 1.0])


        # a dictionary to store the axes objects
        self.ax = {}

        # for displaying the 2D image
        labelkw = dict(fontsize = 5)
        self.ax['2d'] = self.iplot.subplot(1,0)
        plt.setp(self.ax['2d'].get_xticklabels(), **labelkw)
        plt.setp(self.ax['2d'].get_yticklabels(), **labelkw)
        # for displaying cuts along the dispersion direction
        self.ax['slicew'] = self.iplot.subplot(0,0,sharex=self.ax['2d'])
        self.ax['slicew'].set_title(self.aperture.name, fontsize=8)
        plt.setp(self.ax['slicew'].get_xticklabels(), visible=False)
        plt.setp(self.ax['slicew'].get_yticklabels(), **labelkw)
        # for display cuts along the cross-dispersion direction
        self.ax['slices'] = self.iplot.subplot(1,1,sharey=self.ax['2d'])
        self.ax['slices'].xaxis.tick_top()
        plt.setp(self.ax['slices'].get_xticklabels(), rotation=270, **labelkw)
        plt.setp(self.ax['slices'].get_yticklabels(), visible=False)

        # set up some plotting defaults
        self.percentiles = percentiles
        self.extent = [ self.waxis.min(), self.waxis.max(),
                        self.saxis.min(), self.saxis.max()]


        # a dictionary to store plotted objects, so they can be modified
        self.plotted = {}

        # plot the image
        self.plotted['2d'] = self.ax['2d'].imshow(self.imagetoplot,
                                                  cmap='gray',
                                                  extent=self.extent,
                                                  interpolation='nearest',
                                                  aspect='auto',
                                                  zorder=0,
                                                  origin='lower',
                                                  vmin=self.vmin,
                                                  vmax=self.vmax)
        self.ax['2d'].set_xlim(self.extent[0:2])
        self.ax['2d'].set_ylim(self.extent[2:4])

        # add the trace guess
        self.plotted['traceguess'] =  self.ax['2d'].plot([],[],
                                            color='mediumseagreen',
                                            linestyle='--')[0]

        # add the trace points
        self.plotted['tracepoints'] = self.ax['2d'].plot([],[],
                                            color='seagreen',
                                            linewidth=0,
                                            marker='o')[0]

        # add the sky subtraction mask
        self.plotted['skymask'] =  self.ax['2d'].imshow(self.skymaskimage,
                                            cmap=one2another('deepskyblue', 'deepskyblue', alphabottom=0.0, alphatop=1.0),
                                            extent=self.extent,
                                            interpolation='nearest',
                                            aspect='auto',
                                            zorder=100,
                                            vmin=0.5, vmax=1.5,
                                            origin='lower')

        # add the mask flux mask
        #self.plotted['maskfluxmask'] =  self.ax['2d'].imshow(self.maskfluxmaskimage,
        #                                    cmap=one2another('gold', 'gold', alphabottom=0.0, alphatop=1.0),
        #                                    extent=self.extent,
        #                                    interpolation='nearest',
        #                                    aspect='auto',
        #                                    zorder=100,
        #                                    vmin=0.5, vmax=1.5,
        #                                    origin='lower')

        offsets = (np.array([-1,1])[np.newaxis, :]*self.extractionwidths[:,np.newaxis]).flatten()
        x = self.waxis
        y = self.traceCenter(x)[:,np.newaxis] + offsets[np.newaxis, :]
        self.plotted['extractionedges'] = self.ax['2d'].plot(x, y,
                                            linewidth=1,
                                            alpha=0.25,
                                            color='turquoise',
                                            zorder = 1000)

        # add the sky subtraction mask
        self.plotted['extractionmask'] =  self.ax['2d'].imshow(
                                            self.extractionmaskimage,
                                            cmap=one2another('aquamarine', 'aquamarine', alphabottom=0.0, alphatop=1.0),
                                            extent=self.extent,
                                            interpolation='nearest',
                                            aspect='auto',
                                            zorder=100,
                                            vmin=0.5, vmax=1.5,
                                            origin='lower')

        self.labelWidths()

        # add crosshair
        crosskw = dict(alpha=0.5, color='darkorange', linewidth=1)
        self.plotted['crosss'] = self.ax['2d'].axvline(self.crosshair['w'],
                                                        **crosskw)
        self.plotted['crosssextend'] = self.ax['slicew'].axvline(
                                                        self.crosshair['w'],
                                                        linestyle='--',
                                                        **crosskw)
        self.plotted['crossw'] = self.ax['2d'].axhline(self.crosshair['s'],
                                                        **crosskw)
        self.plotted['crosswextend'] = self.ax['slices'].axhline(
                                                        self.crosshair['s'],
                                                        linestyle='--',
                                                        **crosskw)
        # plot slices
        slicekw = dict(color='darkorange', linewidth=1, alpha=1)
        self.plotted['slices'] = self.ax['slices'].plot(*self.slices,
                                                        **slicekw)[0]
        self.plotted['slicew'] = self.ax['slicew'].plot(*self.slicew,
                                                        **slicekw)[0]

        # plot the recently dragged regions
        self.plotted['dragger'] = self.ax['2d'].plot([],[],
                                                        marker='+',
                                                        alpha=0.25,
                                                        color='royalblue')[0]

        # plot the estimated centroids and widths
        for x in ['up', 'down']:
            self.plotted['fittedwidths{0}'.format(x)]=self.ax['2d'].plot([],[],
                                                linewidth=1,
                                                color='mediumorchid',
                                                alpha=1.0)[0]

        # start off unconverged
        self.notconverged = True
        plt.draw()

    def labelWidths(self):
        self.plotted['widthlabels'] = [
                self.ax['2d'].text(
                        (i+0.5)/float(self.numberofapertures), 0.1,
                        'extract\n{:.1f} pixels\nfrom center'.format(self.extractionwidths[i]),
                        fontsize=8,
                        color='aquamarine',
                        alpha=0.5,
                        ha='center',
                        va='center',
                        transform=self.ax['2d'].transAxes)
                for i in range(self.numberofapertures)]

    def moveCrosshair(self, pressed=None, w=None, s=None):
        '''use new values of w and s to move the crosshair and replot'''

        # pull the values from the mouse click
        if pressed is not None:
            w = pressed.xdata
            s = pressed.ydata

        # update the stored value
        if w is not None:
            self.crosshair['w'] = w
        if s is not None:
            self.crosshair['s'] = s

        # update the position on the 2D plot
        self.plotted['crosss'].set_xdata(self.crosshair['w'])
        self.plotted['crosssextend'].set_xdata(self.crosshair['w'])
        self.plotted['crossw'].set_ydata(self.crosshair['s'])
        self.plotted['crosswextend'].set_ydata(self.crosshair['s'])

        # update the slices in the 1D plots
        self.plotted['slices'].set_data(*self.slices)
        self.ax['slices'].set_xlim(0, self.slices[0].max())
        self.plotted['slicew'].set_data(*self.slicew)
        self.ax['slicew'].set_ylim(0, self.slicew[1].max())

    @property
    def slices(self):
        '''return y, x of a slice along the spatial direction'''
        i = int(np.interp(self.crosshair['w'], self.waxis, np.arange(len(self.waxis))))
        return self.aperture.images['science'][i,:], self.saxis

    @property
    def slicew(self):
        '''return x, y of a slice along the wavelength direction'''
        i = int(np.interp(self.crosshair['s'], self.saxis, np.arange(len(self.saxis))))
        return self.waxis, self.aperture.images['science'][:,i]

    @property
    def imagetoplot(self):
        '''for plotting, the science image'''
        return np.transpose(np.log(self.aperture.images['science']))

    def extractionmask(self, width):
        '''to define those pixels that fall within the default extraction mask,
            and, for those along the edges, their fractional weights'''

        distancefromtrace = np.abs(self.aperture.s - self.traceCenter(self.aperture.w))

        # create an empty mask
        mask = np.zeros_like(self.aperture.s)

        # assign 1.0 to all those pixels that are entirely within the aperture
        definitelyin = distancefromtrace <= np.floor(width)
        mask[definitelyin] = 1.0
        # assign fractional weight to those pixels on the border
        fraction = distancefromtrace - np.floor(width)
        border = (fraction > 0)*(fraction < 1)
        mask[border] = 1.0 - fraction[border]

        return mask

    def skymask(self, width):
        '''to define those pixels that are considered sky'''

        # create a blank mask
        mask = np.zeros_like(self.aperture.images['science'])

        # loop through the sky offsets, and use them to add and subtract
        for d in self.skyoffsets:
            top, bottom, whattodo = d['top'], d['bottom'], d['whattodo']
            # identifiy pixels that been selected
            ok  = self.aperture.s > (self.traceCenter(self.aperture.w) + bottom)
            ok *= self.aperture.s < (self.traceCenter(self.aperture.w) + top)
            # either add or remove them from the sky mask
            mask[ok] = whattodo

        absolutedistancefromtrace = np.abs(self.aperture.s - self.traceCenter(self.aperture.w))
        mask[absolutedistancefromtrace < (width + self.instrument.extractiondefaults['skyGap'])] = 0
        #self.speak('')
        #self.speak('{}'.format(width + self.obs.skyGap))
        # return the populated map
        return mask

    #def maskfluxmask(self, width):
        '''to define those pixels that are considered sky'''

        # create a blank mask
    #    mask = np.zeros_like(self.aperture.images['science'])

        # loop through the sky offsets, and use them to add and subtract
    #    for d in self.maskfluxoffsets:
    #        top, bottom, whattodo = d['top'], d['bottom'], d['whattodo']
            # identifiy pixels that been selected
    #        ok  = self.aperture.s > (self.traceCenter(self.aperture.w) + bottom)
    #        ok *= self.aperture.s < (self.traceCenter(self.aperture.w) + top)
            # either add or remove them from the sky mask
    #        mask[ok] = whattodo

    #    absolutedistancefromtrace = np.abs(self.aperture.s - self.traceCenter(self.aperture.w))
    #    mask[absolutedistancefromtrace < (width + 50)] = 0
        #self.speak('')
        #self.speak('{}'.format(width + self.obs.skyGap))
        # return the populated map
    #    return mask


    @property
    def vmin(self):
        '''for plotting, the minimum value, for the imshow grayscale'''
        ok = np.isfinite(self.imagetoplot)
        return np.percentile(self.imagetoplot[ok], self.percentiles[0])

    @property
    def vmax(self):
        '''for plotting, the maximum value, for the imshow grayscale'''
        ok = np.isfinite(self.imagetoplot)
        return np.percentile(self.imagetoplot[ok], self.percentiles[1])

    @property
    def skymaskimage(self):
        '''for plotting, a masked array of the sky mask'''
        image = np.zeros_like(self.imagetoplot)
        chunksize = image.shape[1]/self.numberofapertures

        for i,width in enumerate(self.extractionwidths):
            mask = self.skymask(width).T

            image[:,int(i*chunksize):] = mask[:,int(i*chunksize):]*(1.0 - 0.2*(i%2))

        return image/image.max()

    #@property
    #def maskfluxmaskimage(self):
    #    '''for plotting, a masked array of the mask flux mask'''
    #    image = np.zeros_like(self.imagetoplot)
    #    chunksize = image.shape[1]/self.numberofapertures

    #    for i,width in enumerate(self.extractionwidths):
    #        mask = self.maskfluxmask(width).T

    #        image[:,int(i*chunksize):] = mask[:,int(i*chunksize):]*(1.0 - 0.2*(i%2))

    #    return image/image.max()

    @property
    def extractionmaskimage(self):
        '''for plotting, a masked array of the extraction mask'''
        image = np.zeros_like(self.imagetoplot)
        chunksize = image.shape[1]/self.numberofapertures

        for i,width in enumerate(self.extractionwidths):
            mask = self.extractionmask(width).T
            image[:,int(i*chunksize):] = mask[:,int(i*chunksize):]*(1.0 - 0.2*(i%2))

        return image/image.max()
