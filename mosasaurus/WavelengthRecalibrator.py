'''
The tools in the wavelength recalibrator allow us to shift and stretch
the wavelength calibrations within an unshifted spectral cube, to make
strong features line up across all stars and all exposures.
'''

from .imports import *

def subtractcontinuum(y):
		'''
		Subtract off a rough continuum, so it's easier to correlate.
		'''

		# figure out an overall slope
		x = np.arange(len(y))
		initialfit = np.polyfit(x, y, 1)
		model = np.polyval(initialfit, x)

		# fit only to the things in the top 30% above that slope
		top = y > np.percentile(y - model, 70)
		fit = np.polyfit(x[top], y[top], 1)
		model = np.polyval(fit, x)

		subtracted = y - model
		return subtracted/np.std(subtracted)
		#return residy - np.mean(y)

class WavelengthRecalibrator(Talker):
	'''
	This object latches onto an unshifted cube, and uses it
	to determine the streches needed to line up all its spectra.
	'''

	def __init__(self, cube, visualize=False):
		'''
		initialize by assigning to a cube
		'''

		Talker.__init__(self)
		# make sure that we're starting with an unshifted cube
		assert(cube.shift == False)

		# store the unshifted cube in here
		self.unshiftedcube = cube

		#
		self.visualize = visualize

		# load the correction, or create it from scratch
		try:
			self.load()
			self.speak('loaded shifts successfully')
		except IOError:
			self.determineStretches()
			self.save()

	@property
	def filename(self):
		'''
		What's the filename associated with the spectral stretch?
		'''
		return os.path.join(self.unshiftedcube.directory, 'spectralstretch.npy')

	# keep - save the shifts to a file
	def save(self):
		'''Save the spectral shifts and stretches to a file.'''

		# what's the filename to save to?
		self.speak('exported the spectral shifts and stretches to {}'.format(self.filename))
		tosave = self.corrections
		np.save(self.filename, tosave)

	# keep - load the shifts from a file
	def load(self):
		'''Load the spectral shifts and stretches from a file.'''

		self.speak('attempting to load previously saved shifts and stretches from {}'.format(self.filename))
		loaded = np.load(self.filename)[()]
		self.corrections = loaded
		self.speak('  loading [corrections] from the saved cube structure')

	def determineStretches(self, plot=True):
		'''
		Starting from an unshifted cube, determine a shift and a stretch
		that would be required for all of the spectra, across all wavelengths
		to line up on fixed wavelengths.

		This relies on the existences of strong features in the spectra,
		in wavelength ranges specified in the instrument for the given disperser,
		that are common across all stars.

		SHORT TERM -- make it so that each star could ignore some of the regions,
		if that region either falls outside of its wavelength range or if
		that particular line doesn't exist in the star.

		LONG TERM -- rely on stellar and telluric models for doing this
		cross-matching, so we can be sure to land on absolute wavelength scales.
		This will require quite a bit of development and potentially some
		reorganization.
		'''

		self.speak('shifting all spectra to line up...')
		self.speak('')


		self.corrections = {}
		# store an array of file prefixes associated with each timepoint
		self.corrections['prefixes'] = self.unshiftedcube.obs.fileprefixes['science']
		# what's the wavelength zero-point for the polynomial correction?
		self.corrections['midpoint'] = np.mean(self.unshiftedcube.spectral['wavelength'])
		self.corrections['commonwavelength'] = self.unshiftedcube.spectral['wavelength']


		###subset['midpoint'] = np.mean(self.unshiftedcube.spectral['wavelength'])
		for line in self.unshiftedcube.obs.instrument.alignmentranges:
				linekey = 'offset_{}'.format(line)
				self.corrections[linekey] = {}
				for star in self.unshiftedcube.stars:
						self.corrections[linekey][star] = {}

		for key in ['shift', 'stretch']:
				self.corrections[key] = {}
				for star in self.unshiftedcube.stars:
						self.corrections[key][star] = {}

		# pull out the cube of counts
		c = self.unshiftedcube.cubes['raw_counts']

		# pull out an array of ranges to consider for correlations
		alignmentranges = self.unshiftedcube.obs.instrument.alignmentranges

		# a dictionary to store each line's individual shift
		localshifts = {}

		# create an a structure to store the shifts in
		for star in self.unshiftedcube.stars:
				localshifts[star] = {}
				#shift = np.zeros(self.unshiftedcube.numberoftimes)

		if self.visualize:
			# make a plot to put everything in
			plt.figure(figsize=(20,10), dpi=72)
			gs = plt.matplotlib.gridspec.GridSpec(3*len(self.unshiftedcube.stars) + 1,len(alignmentranges),
									hspace=0.25,
									left=0.05, right=0.95, bottom=0.05, top=0.95, height_ratios=[1,1,0.5]*len(self.unshiftedcube.stars)+[1.5])
			self.axshifts = {}

		# loop over timepoints
		for i, prefix in enumerate(self.corrections['prefixes']):
				if self.visualize:
					for a in self.axshifts.values():
						a.cla()

				for iline, line in enumerate(alignmentranges.keys()):

						# set up the arrays for this correlation
						l, r = alignmentranges[line]

						# define the range for this correlation
						left = np.argmin(np.abs(self.unshiftedcube.spectral['wavelength'] - l))
						right = np.argmin(np.abs(self.unshiftedcube.spectral['wavelength'] - r))

						# create a master template
						wave = self.unshiftedcube.spectral['wavelength'][left:right]
						masterexposure = 0
						master = c[self.unshiftedcube.target][masterexposure, left:right]
						start = subtractcontinuum(master)

						# calculate all the shifts
						# set up the plotting axes

						for istar, star in enumerate(self.unshiftedcube.stars):



							# pull out the spectrum for one star,
							# at one extraction width, at one time point,
							# and trim it to a narrow range around the correlation anchors
							spectrum = c[star][i,left:right]

							# subtract its continuum
							this = subtractcontinuum(spectrum)

							# cross correlation with the anchors
							xc = np.correlate(this, start, 'same')

							# do a quadratic fit to estimate the peak
							x = np.arange(len(xc))
							coeff = np.polynomial.polynomial.polyfit(x[xc.argmax()-5:xc.argmax()+5], xc[xc.argmax()-5:xc.argmax()+5], 2)
							fit = np.polynomial.polynomial.Polynomial(coeff)
							der = fit.deriv()
							peak = der.roots()

							# calculate the offset from the peak
							offset = peak - len(wave)/2 + 1
							linecenter = (l + r)/2.0
							localshifts[star][linecenter] = offset

							# store this for later
							# self.unshiftedcube.squares['offset_{}'.format(line)][star][i] = offset
							self.corrections['offset_{}'.format(line)][star][prefix] = offset

							if self.visualize:
								# plot spectrum and reference
								pair = '{line}+{star}'.format(**locals())
								try:
										ax = self.axshifts['spectrum+'+pair]
								except KeyError:
										ax = plt.subplot(gs[3*istar + 1,iline])
										self.axshifts['spectrum+'+pair] = ax

								ax.set_xlim(wave.min(), wave.max())
								ax.set_autoscaley_on(True)

								ax.plot(wave, start/np.std(start), color='gray', linewidth=3, alpha=0.4)
								ax.plot(wave, this/np.std(this), color=self.unshiftedcube.starcolor(star))
								if iline == 0:
									ax.set_ylabel('{}/{}'.format(i+1, self.unshiftedcube.numberoftimes))
								plt.setp(ax.get_yticklabels(), visible=False)
								ax.set_xlabel('Wavelength\n(angstroms)')


								# plot the correlation function
								try:
									ax = self.axshifts['correlation+'+pair]
								except KeyError:
									ax = plt.subplot(gs[3*istar,iline])
									self.axshifts['correlation+'+pair] = ax
								ax.set_autoscaley_on(True)
								ax.plot(x - len(wave)/2 + 1, xc, alpha=0.5, color='gray')
								plt.setp(ax.get_yticklabels(), visible=False)
								ax.scatter(peak - len(wave)/2 + 1, fit(peak), color='gray')
								#ax.set_xlim(0,len(x))
								if iline == 0:
										ax.set_ylabel('{}'.format(star))

								if istar == 0:
										ax.set_title(line)

								ax.axvline(len(x)/2.0 - len(wave)/2, color='gray', zorder=-2, alpha=0.4, linestyle='--')

				if self.visualize:
					try:
						ax = self.axshifts['compilation']
					except KeyError:
						ax = plt.subplot(gs[-1,:])
						self.axshifts['compilation'] = ax
					plt.sca(ax)
				for star in localshifts.keys():

						# pull out the approximate line centers
						linecenters = list(localshifts[star].keys())

						# pull out the offsets associated with each
						offsets = [localshifts[star][l][0] for l in linecenters]

						# fit a function for offset vs. original wavelength
						originalwavelength = self.unshiftedcube.spectral['wavelength']
						midpoint = np.mean(originalwavelength)
						x = linecenters - midpoint
						fit = np.polyfit(x, offsets, 1)

						# store the corrections (stretch + shift)
						self.corrections['stretch'][star][prefix] = fit[0]
						self.corrections['shift'][star][prefix] = fit[1]


						dw = np.polyval(fit, originalwavelength - midpoint)
						phrase = 'dw = {:.4}xw{:+.4}'.format(*fit)

						self.speak('stretch for {}+{} is {}'.format(star, prefix, phrase))
						# (MOVE THIS TO EXTRACTED TO SUPERSAMPLED!)
						#for key in self.unshiftedcube.cubes.keys():
						#		interpolated = fluxconservingresample(
						#												originalwavelength - dw,
						#												self.unshiftedcube.cubes[key][star][i,:],
						#												originalwavelength)
						#		self.unshiftedcube.cubes[key][star][i,:] = interpolated
						#		self.speak('shifted [{}] by {}'.format(key, phrase))
						if self.visualize:
							color = self.unshiftedcube.starcolor(star)
							plt.plot(self.unshiftedcube.spectral['wavelength'], dw, color=color)
							plt.scatter(linecenters, offsets, color=color)

				if self.visualize:
					plt.xlim(np.min(self.unshiftedcube.spectral['wavelength']), np.max(self.unshiftedcube.spectral['wavelength']))
					plt.xlabel('Original Wavelength (angstroms)')
					plt.ylim(-15, 15)
					plt.ylabel('New - Original')

					pltdir = os.path.join(self.unshiftedcube.directory, 'shifts')
					mkdir(pltdir)
					pltfilename = os.path.join(pltdir, 'shift_{}.pdf'.format(self.unshiftedcube.obs.fileprefixes['science'][i]))
					plt.savefig(pltfilename)

					#self.speak( "shift = {4} for star {0}; {2}/{3} spectra".format(star, self.unshiftedcube.numberoftimes, i+1, len(c[star][:,0]), offset))
					self.speak('plot saved to {}'.format(pltfilename))

	def recreateSupersampled(self):
		'''
		Loop through all the extracted spectra,
		and resample them onto their new grids,
		using the new stretched wavelength calibration
		necessary for each.
		'''

		# these are things where we care about the sum matching extracted to supersampled
		self.additivekeys = ['raw_counts', 'sky']
		# these are things where we want the individual values matching extracted to supersampled
		self.intrinsickeys = ['centroid', 'width', 'peak']
		# these are all the keys that will be supersampled
		self.keys = self.additivekeys + self.intrinsickeys

		# what's the wavelength grid we're resampling onto
		self.commonwavelength = self.unshiftedcube.spectral['wavelength']

		# loop over stars
		for star in self.unshiftedcube.stars:
			directory = os.path.join(self.unshiftedcube.directory, star)

			# loop over prefixes
			for prefix in self.corrections['prefixes']:
				extractedFilename = os.path.join(directory, 'extracted_{}.npy'.format(prefix))

				# load the original extracted file
				extracted = np.load(extractedFilename)

				# nudge the wavelengths
				originalwavelength = extracted['wavelength'] + 0.0
				midpoint = self.corrections['midpoint']
				coefficients = self.corrections['stretch'][star][prefix], self.corrections['shift'][star][prefix]
				dw = np.polyval(coefficients, originalwavelength - midpoint)
				wavelength = originalwavelength + dw
				phrase = 'dw = {:.4}x(w - {midpoint}){:+.4}'.format(*coefficients, midpoint=midpoint)
				self.speak('nudge wavelengths for {} by {}'.format(prefix, phrase))
				extracted['wavelength'] = wavelength

				# save the new shifted spectrum out to a file
				shiftedExtractedFilename = os.path.join(directory, 'shiftedextracted_{}.npy'.format(prefix))
				np.save(shiftedExtractedFilename, extracted)

				# this function should create a supersampled dictionary of a spectrum
				supersampled = extracted2supersampled(
									extracted=extracted,
									wavelength=self.commonwavelength,
									additivekeys=['raw_counts', 'sky'],
									intrinsickeys=['centroid', 'width', 'peak'],
									)
