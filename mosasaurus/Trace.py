'''Trace defines the extraction and sky subtraction regions for an Aperture'''

from imports import *


class Trace(Talker):
    def __init__(self, aperture):
        Talker.__init__(self)
        self.aperture = aperture
        self.obs = self.aperture.obs

        self.saxis = self.aperture.saxis
        self.waxis = self.aperture.waxis

        # pull out values from the aperture
        self.order = self.obs.traceOrder
        self.traceWidth = self.obs.extractionWidth

        self.crosshair = dict(s=0.0, w=0.0)
        o = self.obs
        inner = o.extractionWidth + o.skyGap
        outer = inner + o.skyWidth
        self.skyoffsets = [ dict(top=outer, bottom=inner, whattodo=1),
                            dict(top=-inner, bottom=-outer, whattodo=1)]

        self.tracepoints = {'w':[], 's':[]}
        self.traceguess = np.poly1d(0.0)
        self.setup()
        #self.updateMasks()
        #for i in range(3):
        #    self.fitTrace()
        self.run()
        self.save()

    def save(self):
        filename = self.aperture.directory + 'trace_{0}.npy'.format(self.aperture.name)
        np.save(filename, (self.tracefitcoef, self.tracefitwidth))
        self.speak("saved trace parameters to {0}".format(filename))

        filename = self.aperture.directory + 'trace_{0}.pdf'.format(self.aperture.name)
        self.figure.savefig(filename)

        filename = self.aperture.directory + 'skyMask_{0}.npy'.format(self.aperture.name)
        np.save(filename, self.skymask)
        self.speak("saved a sky mask to {0}".format(filename))

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

                    press and hold "e" at a point on the sky,
                    and release at a different point on the sky

                        "Extend" the sky subtraction region.
                        All this cares about is the offset
                        in the cross-dispersion direction
                        away from the spectral trace.

                    press and hold "r" at a point on the sky,
                    and release at a different point on the sky

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
        stillgoing = True
        while stillgoing:
            self.speak('Please refine the extraction aperture.')
            self.speak(' Your options are:')
            self.speak('  [q] save and quit')
            self.speak('  [c] move the crosshair, and plot slices along it')
            self.speak('  [t] add a trace point, somewhere along the star')
            self.speak('  [e] extend a sky region (twice to define the start and stop)')
            self.speak('  [r] remove area from a sky region (twice to start and stop)')
            self.speak("  [f] to fit the trace, using the star's centroids and the sky areas")
            pressed = self.iplot.getKeyboard()
            if pressed.key == 'q':
                stillgoing = False
                break
            elif pressed.inaxes == None:
                continue
            elif pressed.key == 'c':
                self.moveCrosshair(pressed)
            elif pressed.key == 't':
                self.addTracePoint(pressed)
                self.moveCrosshair(pressed)
            elif pressed.key == 'e':
                self.modifySky(pressed)
            elif pressed.key == 'r':
                self.modifySky(pressed)
            elif pressed.key == 'f':
                self.fitTrace()
            else:
                self.speak('Hmmm, nothing has be defined for "{}"...'.format(pressed.key))
            plt.draw()

    @property
    def traceCenter(self):
        try:
            return self.tracefit
        except AttributeError:
            return self.traceguess

    def modifySky(self, pressed, whattodo=None):
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
            print d
            try:
                self.skyoffsets.append(d)
            except AttributeError:
                self.skyoffsets = [d]

            self.updateMasks()
            #ma = np.ma.MaskedArray(self.images['Science'], mask==0)
            #skyspectrum = np.ma.median(ma, self.aperture.sindex)

            # remove the dictionary, to start over again for the next
            del self.skykeys
        elif len(self.skykeys) == 1:
            self.speak('Please do another "{}" again to edit the sky.'.format(pressed.key))
        else:
            self.speak('How did you get here?')

    def updateMasks(self):
        '''update the sky mask (because either the sky offsets have changed),
            or the trace center has changed), and update the plots'''

        self.plotted['extractionmask'].set_data(self.extractionmaskimage)
        self.plotted['skymask'].set_data(self.skymaskimage)

        plt.draw()

    def addTracePoint(self, pressed):
        '''from a KeyEvent, add a point to the list of trace guesses'''
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
        print xfine
        print self.traceguess(xfine)

        # the sky moves, if the trace does
        self.updateMasks()

    @property
    def images(self):
        return self.aperture.images

    def fitTrace(self):
        '''using the guess of the trace, and the sky subtraction regions,
            fit for trace using actual centroids'''



        # estimate a rough 1D spectrum
        flattened = self.images['Science']/self.images['NormalizedFlat']
        roughSky1d=np.average(flattened,
                                axis=self.aperture.sindex,
                                weights=self.skymask)
        reshaped = roughSky1d.reshape((self.waxis.shape[0],1))
        self.images['Sky'] = reshaped*np.ones_like(flattened)
        self.images['Subtracted'] = flattened - self.images['Sky']
        fine = np.isfinite(self.images['Subtracted'])

        considerstar = self.extractionmask

        weights = np.maximum(fine*considerstar*self.images['Subtracted'], 0)
        weights += 1.0/np.sum(weights)

        # use the roughly subtracted image to fit centroids
        flux = np.average(self.images['Subtracted'],
                            axis=self.aperture.sindex,
                            weights=weights)

        if (np.sum(weights, self.aperture.sindex) <= 0).any():
            self.speak('weights were wonky')
            assert(False)

        try:
            fluxWeightedCentroids = np.average(self.aperture.s,
                        axis=self.aperture.sindex,
                        weights=fine*considerstar*self.images['Subtracted'])
        except ZeroDivisionError:
            self.speak('UH-OH, got zero-division error')
            return
        if np.isfinite(fluxWeightedCentroids).all() == False:
            self.speak('centroids were wacky')
            assert(False)
        reshapedFWC = fluxWeightedCentroids[:,np.newaxis]
        weights = fine*considerstar*self.images['Subtracted']
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
        #self.traceWidth = width
        self.updateMasks()

        ''' # create a rough LSF
        self.images['RoughLSF'] = np.exp(-0.5*((self.s - self.traceCenter(self.w))/self.traceWidth)**2)

        converged = np.abs(fit_width - old_width) < 0.01
        self.speak( "  {0} -> {1}, converged = {2}".format(old_width, width, converged))
        old_width = width

        #self.input('continue?')
        attempts += 1

        if attempts > 10:
        converged = True

        self.traceCenter = np.polynomial.polynomial.Polynomial(self.tracefitcoef)
        self.traceWidth = width'''



    def setup(self, percentiles=[10,90]):
        '''make the plotting window and interactive tools'''

        plt.ion()
        self.figure = plt.figure('tracing the spectrum',
                                figsize=(8,4), dpi=100)

        # create an interactive plot
        self.iplot = zachopy.iplot.iplot(2,2,
                                        hspace=0, wspace=0,
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
        self.plotted['2d'] = self.ax['2d'].imshow(self.image,
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
                                            alpha=0.5,
                                            cmap='winter_r',
                                            extent=self.extent,
                                            interpolation='nearest',
                                            aspect='auto',
                                            zorder=100,
                                            vmin=0.5, vmax=1.5,
                                            origin='lower')

        # add the sky subtraction mask
        self.plotted['extractionmask'] =  self.ax['2d'].imshow(
                                            self.extractionmaskimage,
                                            alpha=0.5,
                                            cmap='summer_r',
                                            extent=self.extent,
                                            interpolation='nearest',
                                            aspect='auto',
                                            zorder=100,
                                            vmin=0.5, vmax=1.5,
                                            origin='lower')

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
                                                color='darkcyan',
                                                alpha=1.0)[0]

        plt.draw()

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
        i = np.interp(self.crosshair['w'], self.waxis, np.arange(len(self.waxis)))
        return self.images['Science'][i,:], self.saxis

    @property
    def slicew(self):
        '''return x, y of a slice along the wavelength direction'''
        i = np.interp(self.crosshair['s'], self.saxis, np.arange(len(self.saxis)))
        return self.waxis, self.images['Science'][:,i]

    @property
    def image(self):
        return np.transpose(np.log(self.images['Science']))

    @property
    def vmin(self):
        ok = np.isfinite(self.image)
        return np.percentile(self.image[ok], self.percentiles[0])

    @property
    def vmax(self):
        ok = np.isfinite(self.image)
        return np.percentile(self.image[ok], self.percentiles[1])

    @property
    def extractionmask(self):
        fromtrace = self.aperture.s - self.traceCenter(self.aperture.w)
        return np.abs(fromtrace) < self.traceWidth

    @property
    def skymask(self):
        mask = np.zeros_like(self.images['Science'])
        for d in self.skyoffsets:
            top, bottom, whattodo = d['top'], d['bottom'], d['whattodo']
            # modify the sky mask
            ok  = self.aperture.s > (self.traceCenter(self.aperture.w) + bottom)
            ok *= self.aperture.s < (self.traceCenter(self.aperture.w) + top)
            mask[ok] = whattodo
        return mask

    @property
    def skymaskimage(self):
        skymasktoplot = np.ma.masked_where(self.skymask == 0, self.skymask)
        return np.transpose(skymasktoplot)

    @property
    def extractionmaskimage(self):
        extractionmasktoplot = np.ma.masked_where(self.extractionmask == 0, self.extractionmask)
        return np.transpose(extractionmasktoplot)

    """
              # try refitting with a better initial guess!
              self.createTrace(inputCoef=coef)


          def displayTrace(self):
              self.display.rgb( self.images['Subtracted'],
                                self.images['RoughLSF'],
                                self.considersky)


          def createSkyApertures(self, visualize=True):
            '''Let user select the best sky apertures.'''
            self.speak("setting up the sky apertures")
            filename = self.directory + 'skyMask_{0}.npy'.format(self.aperture.name)
            try:
              self.images['skyMask'] = np.load(filename)
              self.speak("loaded sky apertures from {0}".format(filename))
            except IOError:
              finished = False
              plt.figure('sky apertures', figsize=(10,10), dpi=50)
              i = zachopy.iplot.iplot(2,2,
                                        hspace=0, wspace=0,
                                        height_ratios=[.2, 1], width_ratios=[1, .2])

              self.aximage = i.subplot(1,0)
              self.axskyspectrum = i.subplot(0,0, sharex=self.aximage)
              self.axskyprofile = i.subplot(1,1, sharey=self.aximage)
              axes = [self.aximage, self.axskyspectrum, self.axskyprofile]
              for ax in axes:
                  plt.setp(ax.get_xticklabels(), visible=False)
                  plt.setp(ax.get_yticklabels(), visible=False)

              mask = np.zeros_like(self.images['Science'])
              extent=[self.waxis.min(), self.waxis.max(), self.saxis.min(), self.saxis.max()]
              finishedplotting = True
              first = True
              while(finished == False):


                if first == False:
                  # have user select a sky region
                  self.speak("please click to select a sky region")
                  clicks = i.getMouseClicks(n=2)

                  # clear the axes
                  for a in axes:
                    a.cla()


                # display the image
                values = np.percentile(self.images['Science'], [10,90])
                self.aximage.imshow(np.transpose(np.log(self.images['Science'])),
                              cmap='gray', \
                              extent=extent, \
                              interpolation='nearest', aspect='auto', \
                              vmin=np.log(values[0]), vmax=np.log(values[1]))
                self.aximage.imshow(np.transpose(mask), alpha=0.1, cmap='winter_r', \
                            extent=extent, \
                            interpolation='nearest', aspect='auto')

                # overlay the trace
                self.aximage.plot(self.waxis, self.traceCenter(self.waxis), color='blue', alpha=0.3, linewidth=4)
                self.aximage.set_xlim(self.waxis.min(), self.waxis.max())


                if first == False:
                  # calculate offsets from the trace
                  offsets = (clicks[0].ydata - self.traceCenter(clicks[1].xdata), clicks[1].ydata - self.traceCenter(clicks[1].xdata))
                  bottom = np.min(offsets)
                  top = np.max(offsets)

                  # display the most recent
                  self.aximage.plot(self.waxis, self.traceCenter(self.waxis) + bottom , color='green', alpha=0.3, linewidth=4)
                  self.aximage.plot(self.waxis, self.traceCenter(self.waxis) + top , color='green', alpha=0.3, linewidth=4)
                  mask[(self.s > self.traceCenter(self.w) + bottom) * (self.s < self.traceCenter(self.w) + top)] += 1.0
                  mask = mask > 0

                  ma = np.ma.MaskedArray(self.images['Science'], mask==0)
                  skyspectrum = np.ma.median(ma, self.aperture.sindex)
                  self.axskyspectrum.plot(skyspectrum)
                  click = clicks[-1]
                  self.axskyprofile.cla()
                  self.axskyprofile.plot(self.images['Science'][click.xdata,:], self.saxis)
                  self.axskyprofile.plot((self.images['Science']*mask)[click.xdata,:], self.saxis, linewidth=3)
                  self.axskyprofile.set_xlim((self.images['Science']*mask)[click.xdata,:].min(), (self.images['Science']*mask)[click.xdata,:].max()*2)
                  self.axskyprofile.set_ylim(self.saxis.min(), self.saxis.max())
                  plt.draw()
                  self.speak("Are you happy with the sky subtraction apertures? (default = no)")
                  answer = self.input("  (y)es, (n)o, (r)edo")
                  if "y" in answer:
                    finished = True
                  elif "r" in answer:
                    mask *= 0
                  else:
                    finished = False
                first = False
              self.images['skyMask'] = mask
              np.save(filename, self.images['skyMask'])
              self.speak("saved a sky mask to {0}".format(filename))"""
