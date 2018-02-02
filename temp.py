  def interpolate(self, remake=False, shift=False):
        '''Interpolate the spectra onto a common (uniform) wavelength scale.'''


        # these are things where we care about the sum matching extracted to supersampled
        self.additivekeys = ['raw_counts', 'sky']
        # these are things where we want the individual values matching extracted to supersampled
        self.intrinsickeys = ['centroid', 'width', 'peak']
        # these are all the keys that will be supersampled
        self.keys = self.additivekeys + self.intrinsickeys

        try:
            assert(remake == False)
            self.supersampled = np.load(self.supersampledFilename)
            self.speak('loaded supersampled spectrum from {0}'.format(self.supersampledFilename))
        except (AssertionError, IOError):
            self.speak('creating a supersampled spectrum for {}'.format(self.supersampledFilename))


            self.supersampled = extracted2supersampled(extracted, newwavelength, additivekeys, intrinsickeys)


            # define an empty axes that we're going to populate with spectra
            if self.visualize:
                suptitletext = '{}, {}'.format(self.name, self.exposureprefix)
                try:
                    self.suptitle_supersampled.set_text(suptitletext)
                    self.plotted_supersampled
                except AttributeError:
                    # set up a plot window to show how the interpolation is going
                    self.figure_supersampled = plt.figure('interpolating spectra', figsize=(20,10), dpi=60)

                    self.ax_supersampled = {}
                    sharex=None
                    self.suptitle_supersampled = plt.suptitle(suptitletext)

                    # make subplot for each key, and each aperture
                    gs = plt.matplotlib.gridspec.GridSpec(nkeys, napertures, hspace=0.1, wspace=0.05)
                    self.plotted_supersampled = {}

                try:
                    self.limits_supersampled
                except AttributeError:
                    # for plotting, we'll want to keep track of limits for each key
                    self.limits_supersampled = {}


            # set up a fine, common wavelength grid onto which everything will be interpolated
            try:
                self.supersampled['wavelength']
            except (AttributeError, KeyError):



            # loop over the measurement types
            sharex=None
            for i, key in enumerate(self.keys):

                # supersample onto the grid
                for j,width in enumerate(widths):

                    # make a combined key that includes both the key name, and the widthkey
                    widthkey = '{:04.1f}px'.format(width)
                    combinedkey = key + '_' + widthkey

                    # use flux-conserving resampling to put onto this grid
                    yoriginal = self.extracted[width][key]
                    isnanoriginal = np.isnan(yoriginal)

                    # this ones where the supersampled array was corrupted by nans
                    yclosetonan = fluxconservingresample(
                                        wavelength,
                                        isnanoriginal,
                                        self.supersampled['wavelength']) > 0

                    # supersample onto the uniform grid
                    ysupersampled = fluxconservingresample(
                                        wavelength,
                                        self.extracted[width][key],
                                        self.supersampled['wavelength'],
                                        treatnanas=0.0)

                    assert(np.isfinite(ysupersampled).all())

                    # turn the bad elements back to nans
                    #ysupersampled[yclosetonan] = np.nan

                    if key in self.additivekeys:
                        self.supersampled[combinedkey] = ysupersampled
                    elif key in self.intrinsickeys:
                        self.supersampled[combinedkey] = ysupersampled/self.supersampled['fractionofapixel']
                    else:
                        self.speak("Yikes! It's not clear if {} is an additive or intrinsic quantity!".format(key))

                    # set up the plots to visualize this
                    if self.visualize:
                        try:
                            # does this subplot already exist and just need to be cleared?
                            self.ax_supersampled[combinedkey]
                        except KeyError:
                            # or do we need to create it from scratch?
                            self.ax_supersampled[combinedkey] = plt.subplot(gs[i,j], sharex=sharex)
                            sharex = self.ax_supersampled[combinedkey]

                    # plot demonstration
                    if self.visualize:
                        if key in self.additivekeys:
                            ontheoriginalpixels = self.supersampled[combinedkey]/self.supersampled['fractionofapixel']
                        elif key in self.intrinsickeys:
                            ontheoriginalpixels = self.supersampled[combinedkey]
                        else:
                            ontheoriginalpixels = None

                        try:
                            # can we just updated existing plots?
                            plot_extracted, plot_supersampled = self.plotted_supersampled[combinedkey]
                            plot_extracted[0].set_data(wavelength, self.extracted[width][key])
                            plot_supersampled[0].set_data(self.supersampled['wavelength'], ontheoriginalpixels)
                        except KeyError:
                            # or do we have to make new ones?

                            # plot the original spectrum
                            if key in self.additivekeys:
                                ylabel = combinedkey + "\n(per original pixel)"
                            else:
                                ylabel = combinedkey
                            self.ax_supersampled[combinedkey].set_ylabel(ylabel)
                            plot_extracted = self.ax_supersampled[combinedkey].plot(wavelength, self.extracted[width][key], color='black', alpha=0.5)
                            plot_supersampled = self.ax_supersampled[combinedkey].plot(self.supersampled['wavelength'], ontheoriginalpixels, color='red', alpha=0.5)
                            self.plotted_supersampled[combinedkey] = [plot_extracted, plot_supersampled]
                            self.ax_supersampled[combinedkey].set_xlim(np.min(self.supersampled['wavelength'])-200, np.max(self.supersampled['wavelength']) + 200)
                            try:
                                self.limits_supersampled[combinedkey]
                            except:
                                ok = np.isfinite(ontheoriginalpixels)

                                lims = np.percentile(ontheoriginalpixels[ok], [5, 95])#np.nanmin(), np.nanmax(self.extracted[width][key])
                                span = np.abs(lims[1] - lims[0])
                                nudge = 0.5
                                self.limits_supersampled[combinedkey] = [(lims[0] - span*nudge),  (lims[1] + span*nudge)]
                                if 'centroid' not in key:
                                    self.limits_supersampled[combinedkey][0] = np.maximum(0, self.limits_supersampled[combinedkey][0])
                            self.ax_supersampled[combinedkey].set_ylim(*self.limits_supersampled[combinedkey])

                            if key == self.keys[0]:
                                plt.title(widthkey)
                            if key == self.keys[-1]:
                                plt.xlabel('Wavelength (angstroms)')
                            else:
                                plt.setp(self.ax_supersampled[combinedkey].get_xticklabels(), visible=False)
                            if width != widths[0]:
                                plt.setp(self.ax_supersampled[combinedkey].get_yticklabels(), visible=False)
                                self.ax_supersampled[combinedkey].set_ylabel('')

            if self.visualize:
                #plt.show()
                #plt.draw()
                supersamplingPlotFilename = self.supersampledFilename.replace('npy', 'pdf')
                plt.savefig(supersamplingPlotFilename)
                self.speak('saved a plot to {}'.format(supersamplingPlotFilename))
                #a = self.input('is this okay?')

            #self.input('do you like the interpolation for {0}?'.format(self.name))
            np.save(self.supersampledFilename, self.supersampled)
            self.speak('saved supersampled spectrum to {0}'.format(self.supersampledFilename))
        return self.supersampled

def extracted2supersampled(extracted, newwavelength, additivekeys, intrinsickeys):
    '''
    This function takes a raw extracted spectrum (in pixel coordinates)
    and supersamples it onto another wavelength grid.

    Additive keys are ones that should have the same sum over a given
    wavelength range. Intrinsic keys are ones that should have the same
    values, even if that means the sum would change.
    '''

    # what are all the keys to deal with?
    keys = additivekeys + intrinsickeys

    # we're going to supersample multiple keys, to keep everything together
    nkeys = len(keys)

    # what are the extraction widths to consider?
    widths = np.sort([k for k in extracted.keys() if type(k) != str])
    napertures = len(widths)

    # pull out the wavelength
    wavelength = extracted['wavelength']
    pixelnumber = extracted['w']

    # pull out a reasonable common wavelength grid for this instrument
    commonwavelength = newwavelength

    # what's the (assumed constant) dw for the uniform grid?
    scale = 1.0
    assert((np.diff(commonwavelength) == scale).all())

    # what fraction of an original pixel went into this new pixel?
    doriginaldnew = fluxconservingresample(wavelength, np.ones_like(pixelnumber), commonwavelength)

    # create a supersampled dictionary
    supersampled = {}
    supersampled['wavelength'] = commonwavelength
    supersampled['fractionofapixel'] = doriginaldnew

    # loop over all the keys
    for i, key in enumerate(self.keys):

        #  loop over all the widths
        for j,width in enumerate(widths):

            # make a combined key that includes both the key name, and the widthkey
            widthkey = '{:04.1f}px'.format(width)
            combinedkey = key + '_' + widthkey

            # use flux-conserving resampling to put onto this grid
            yoriginal = extracted[width][key]
            isnanoriginal = np.isnan(yoriginal)

            # this ones where the supersampled array was corrupted by nans
            yclosetonan = fluxconservingresample(
                                wavelength,
                                isnanoriginal,
                                supersampled['wavelength']) > 0

            # supersample onto the uniform grid
            ysupersampled = fluxconservingresample(
                                wavelength,
                                extracted[width][key],
                                supersampled['wavelength'],
                                treatnanas=0.0)

            # make sure the resampling worked OK
            assert(np.isfinite(ysupersampled).all())

            # turn the bad elements back to nans
            ysupersampled[yclosetonan] = np.nan

            if key in additivekeys:
                supersampled[combinedkey] = ysupersampled
            elif key in intrinsickeys:
                supersampled[combinedkey] = ysupersampled/supersampled['fractionofapixel']

    return supersampled
