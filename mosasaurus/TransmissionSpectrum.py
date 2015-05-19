from Observation import Observation
from WavelengthBin import WavelengthBin
from Cube import Cube
from transit import *
import limbdarkening.phoenix as LD
from imports import *
import multiprocessing

class TransmissionSpectrum(Talker):
	""" Transmission spectrum object, which can contain depths + uncertainties, lightcurves, covariance matrices, and structures for refitting every point."""

	def __init__(self, obs=None, binsize=100, remake=False, label='fixedGeometry', maskname='defaultMask'):
		'''Initialize a transmission spectrum object, using a normal Observation file.'''
		Talker.__init__(self, line=200)

		self.speak('initializing a transmission spectrum')

		# should we remake everything?
		self.remake = remake

		# manage the directories
		self.initializeFromObs(obs)

		# manage the bins
		self.setBinsize(binsize)

		# load the light curves associate with this transmission spectrum
		self.constructBins()

		self.speak('its name is {0}'.format(self))

		# manage the mask, and create a default dummy mask if need be
		self.maskname = maskname

		# manage the label, initializing the fit
		self.label = label


	def __repr__(self):
		"""How should this object be represented (e.g. when shown as an element in a list)"""

		try:
			return '<TransmissionSpectrum {name}|{night}|{left}to{right}{unitstring}|{binsize}{unitstring}>'.format(left=self.bins[0].left/self.unit, right=self.bins[-1].right/self.unit, unitstring=self.unitstring, binsize=self.binsize/self.unit, name=self.obs.name,night=self.obs.night)
		except AttributeError:
			return '<TransmissionSpectrum {name}|{night}|[bins unspecified]>'.format(name=self.obs.name,night=self.obs.night)

	def setBinsize(self, binsize, unit=10, unitstring='nm'):
		'''Set the transmission spectrum's binsize.'''
		self.binsize = binsize
		self.unit = unit
		self.unitstring = unitstring

	def initializeFromObs(self, obs):
		'''Initialize a spectrum from an observation filename or instance.'''

		# load the observation structure for this file
		if type(obs) == str:
			self.obs = Observation(obs, nods9=True)
		else:
			self.obs = obs

		# keep track of the name and binsize
		self.name = self.obs.name


	@property
	def binningdirectory(self):
		bd =  self.obs.extractionDirectory + "chromatic{binsize:05.0f}/".format(binsize=self.binsize)
		zachopy.utils.mkdir(bd)
		return bd

	@property
	def rawlightcurvedirectory(self):
		rld = self.binningdirectory + 'originalLCs/'
		return rld

	@property
	def lightcurvedirectory(self):
		ld = self.binningdirectory + 'processedLCs/'
		zachopy.utils.mkdir(ld)
		return ld

	@property
	def fitdirectory(self):
		fd = self.maskdirectory + '{label}/'.format(label=self.label)
		zachopy.utils.mkdir(fd)
		return fd

	@property
	def maskdirectory(self):
		md = self.binningdirectory + "{maskname}/".format(maskname=self.maskname)
		zachopy.utils.mkdir(md)
		return md

	def constructBins(self):
		try:
			# first try loading light curves that have been saved into the directory structure
			possibleTLCs = glob.glob(self.lightcurvedirectory + '*to*/TLC.npy')
			assert(len(possibleTLCs) > 0)
			chunks = []
			for file in possibleTLCs:
				chunks.append(file.split('/')[-2])

		except:
			self.speak("trying to read raw chromatic light curves from {0}".format(self.rawlightcurvedirectory))
			possibleTLCs = glob.glob(self.rawlightcurvedirectory + '*.lightcurve')
			if len(possibleTLCs) == 0:
				self.speak("couldn't find any raw chromatic light curves; trying to recreate them from the spectroscopic cube")
				cube = Cube(self.obs)
				cube.makeLCs(binsize=self.binsize)
				possibleTLCs = glob.glob(self.rawlightcurvedirectory + '*.lightcurve')

			assert(len(possibleTLCs) > 0)
			chunks = []
			for file in possibleTLCs:
				chunks.append((file.split('/')[-1].split('.lightcurve')[0]).split('_')[-1])

		bins = []
		wavelengths = []
		for trimmed in chunks:
			left = np.int(trimmed.split('to')[0])
			right = np.int(trimmed.split('to')[-1])
			try:
				bins.append(WavelengthBin(self, left=left, right=right))
				wavelengths.append(bins[-1].wavelength)
			except IOError:
				pass
		self.bins = np.array(bins)[np.argsort(wavelengths)]
		self.wavelengths = np.array(wavelengths)[np.argsort(wavelengths)]

		# construct a dictionary to link wavelength to list index
		self.w2bin = {}
		for i in range(len(self.wavelengths)):
			self.w2bin[self.wavelengths[i]] = i

		self.bins[0].readTLC()
		# read the first bin (necessary for setting up initial conditions for fits)
		self.speak("spectrum contains {0} bins covering {1} to {2}{3} at {4}{3} resolution".format(len(self.bins), self.bins[0].left/self.unit, self.bins[-1].right/self.unit, self.unitstring, self.binsize/self.unit))


	def loadLCs(self):
		self.speak('loading all the light curves for all the bins''')
		for b in self.bins:
			b.readTLC()

	@property
	def nBins(self):
		"""Return the number of bins in this spectrum."""
		try:
			return self._nBins
		except:
			self._nBins = len(self.bins)
			return self._nBins

	@property
	def nTimes(self):
		"""Return the *maximum* number of times among any of the lightcurves associated with this spectrum."""
		try:
			return self._nTimes
		except:
			# figure out the maximum number of times there are
			self._nTimes = 0
			for i in range(self.nBins):
				b = self.bins[i]
				try:
					b.tlc.flux
				except:
					b.readTLC()
				l = len(b.tlc.flux)
				# update the maximum nTimes, and the most densely populated times (necessary for interpolation)
				if l > self._nTimes:
					self._nTimes = l
			return self._nTimes

	def toArray(self, key):
		'''Create an image array of a TLC variable.'''



		# create an empty array
		a = np.zeros((self.nBins, self.nTimes))

		# we can't guarantee the tlcs all have the exact same timestamps
		#  because some data points might have been thrown out, so we have
		#  to do the interpolation at least once (we can store interpolation
		#  indices inside the bin object to speed this in the future)

		# loop over the bins
		for i in range(self.nBins):

			# pick out this bin's TLC
			tlc = self.bins[i].tlc


			try:
				a[i,:] = tlc.__dict__[key]
			except KeyError:
				a[i,:] = tlc.externalvariables[key]
		return a

	#def interpolationindices(self, tlc):
	#		# define the interpolation indices for this TLC
	#		try:
	#			return tlc.interpolationindices
	#		except:
	#			# interpolate from times to indices
	#			interpolation = scipy.interpolate.interp1d(self.themostdenselypopulatedtimes, np.arange(self.nTimes), bounds_error=True, kind='nearest')
	#			tlc.interpolationindices = interpolation(tlc.bjd).astype(np.int)
	#			return tlc.interpolationindices

	def createMask(self, empty=False, afterfastfit=False):

		if empty:
			maskname = 'defaultMask'
			mask = np.array(self.toArray('bad')).astype(np.byte)
		elif afterfastfit:
			maskname = 'trimOutliers'
			mask = np.array(self.toArray('bad')).astype(np.byte)
		else:
			a = raw_input("What would you like to call this mask?\n ")
			maskname = a.replace(" ","")

			keys = ['peak_', 'sky_', 'width_', 'centroid_']
			# one column for light curves, one for target, one for comparison
			nColumns = 3
			nRows = len(keys)
			plt.figure('masking')
			ip = zachopy.iplot.iplot(nRows, nColumns)

			kw = dict(cmap='gray', interpolation='nearest', aspect='auto', vmin=None, vmax=None)


			mask = np.array(self.toArray('bad')).astype(np.byte)
			nBins = len(self.bins)

			threshold = 4
			# identify outliers
			try:
				for i in range(nBins):
					residuals = self.bins[i].tlc.residuals()
					noise = np.std(residuals)
					bad = np.abs(residuals) > threshold*noise
					print " in bin {0}, {1} points exceeded {2} x {3}".format(i, np.sum(bad), threshold, noise)
					mask[i,:] = mask[i,:] | bad*self.bins[0].tlc.flags['outlier']
			except:
				pass


			# identify "saturated"
			#saturationthreshold = 1.8e6
			#saturated =  (self.toArray('peak_target') > saturationthreshold) |  (self.toArray('peak_target') > saturationthreshold)
			#print "  {0} total points were over {1} in their peak counts".format(np.sum(saturated), saturationthreshold)
			#mask = mask | saturated*self.bins[0].tlc.flags['saturation']


			# set up the axes (this can be slow, so do it once)
			ax = ip.subplot(0,0,name='flux')
			sub = {}
			sub['sharex'], sub['sharey'] = ip.axes['flux'],ip.axes['flux']
			ax = ip.subplot(1,0,name='instrumental',**sub)
			ax = ip.subplot(2,0,name='corrected',**sub)
			ax = ip.subplot(3,0,name='residuals',**sub)
			for i in range(len(keys)):
				key = keys[i]  + 'target'
				ax = ip.subplot(i,1, name=key,**sub)
				key = keys[i]  + 'star01'
				ax = ip.subplot(i,2, name=key,**sub)




			keepgoing = True
			while keepgoing:
				# clear the axes again
				for k in ip.axes.keys():
					ip.axes[k].cla()

				# plot the initial light curve
				flux = self.toArray('flux')
				kw['vmin'], kw['vmax'] = 0.98, 1.02
				ip.axes['flux'].imshow(flux, **kw)
				ip.axes['flux'].set_title('flux')

				# plot the instrumental correction
				instrumental = np.zeros_like(flux)
				residuals = np.zeros_like(flux)
				try:
					for i in range(nBins):
						instrumental[i,:] = self.bins[i].tm.instrument_model()
						residuals[i,:] = self.bins[i].tlc.residuals()
				except:
					pass

				ip.axes['instrumental'].imshow(instrumental, **kw)
				ip.axes['instrumental'].set_title('instrumental')

				ip.axes['corrected'].imshow(flux/instrumental, **kw)
				ip.axes['corrected'].set_title('corrected')

				kw['vmin'], kw['vmax'] = 0.98-1, 1.02-1
				ip.axes['residuals'].imshow(residuals, **kw)
				ip.axes['residuals'].set_title('residuals')


				# plot target diagnostics
				for i in range(len(keys)):
					key = keys[i]  + 'target'
					ax = ip.axes[key]
					ax.set_title(key)
					array = self.toArray(key)
					kw['vmin'], kw['vmax'] =None, None
					ax.imshow(array, **kw)

				# plot comparison diagnostics
				for i in range(len(keys)):
					key = keys[i]  + 'star01'
					ax = ip.axes[key]
					ax.set_title(key)
					array = self.toArray(key)
					kw['vmin'], kw['vmax'] =None, None
					ax.imshow(array, **kw)

				# zoom out a tiny bit to make selection easier
				for k in ip.axes.keys():
					y,x =  self.toArray('flux').shape
					ip.axes[k].set_xlim(0 - x/20, x + x/20)
					ip.axes[k].set_ylim(0 - y/20, y + y/20)


				masked = np.ma.masked_where(mask == 0, mask)
				my_cmap = copy.copy(plt.cm.get_cmap('autumn')) # get a copy of the gray color map
				my_cmap.set_bad(alpha=0)
				for k in ip.axes.keys():
					ax = ip.axes[k]
					ax.imshow(masked, cmap=my_cmap, alpha=0.5, aspect='auto', interpolation='nearest')

				plt.draw()
				unreasonableanswer = True
				while unreasonableanswer:
					answer = raw_input("What would you like to do? [a]dd masking, [s]ubtract masking, [r]efit using this mask, [f]inish?\n   (I'd like to) ")
					unreasonableanswer = False
					if 'a' in answer:
						print "  Click at the two corners of a box you'd like to mask.\n"
						clicks = ip.getMouseClicks(2)

						rows = (clicks[0].ydata, clicks[1].ydata)
						top = np.int(np.max(rows))
						bottom = np.maximum(np.int(np.min(rows)), 0)

						columns =  (clicks[0].xdata, clicks[1].xdata)
						right = np.int(np.max(columns))
						left = np.maximum(np.int(np.min(columns)), 0)


						mask[bottom:top, left:right] = mask[bottom:top, left:right] | self.bins[0].tlc.flags['custom']

					elif 's' in answer:
						print "  Click at the two corners of a box you'd like to unmask.\n"
						clicks = ip.getMouseClicks(2)

						rows = np.round((clicks[0].ydata, clicks[1].ydata))
						top = np.int(np.max(rows))
						bottom = np.maximum(np.int(np.min(rows)), 0)

						columns =  np.round((clicks[0].xdata, clicks[1].xdata))
						right = np.int(np.max(columns))
						left = np.maximum(np.int(np.min(columns)), 0)

						mask[bottom:top, left:right] -= mask[bottom:top, left:right] & self.bins[0].tlc.flags['custom']

					elif 'r' in answer:
						print "Okay, refitting. It may take a while!"
						self.mask = mask
						self.applyMask(maskname=maskname)
						self.fitRigid()
					elif 'f' in answer:
						keepgoing = False
					else:
						unreasonableanswer = True
						print "  I'm sorry, I didn't quite understand that."


		# set the mask and maskname attributes
		self.maskname = maskname
		self.mask = mask

		# save the mask to the mask directory
		self.saveMask()

		# loop through the light curves and apply the mask to each
		self.applyMask(maskname=maskname)

	def saveMask(self):
		"Save the masking array for this particular mask."
		self.speak('saving mask "{0}" to {1}'.format(self.maskname, self.maskdirectory))
		np.save(self.maskdirectory + 'mask.npy', self.mask)

	def applyMask(self, maskname=None):
		"Load a mask from this mask directory and apply it to the light curves."

		# update the maskname attribute to the desired one
		if maskname is not None:
			self.maskname = maskname

		# load the mask from the masking directory (probably slows things down a bit, but ensure everything links up)
		try:
			self.mask = np.load(self.maskdirectory + 'mask.npy')
		except:
			self.speak("couldn't load requested mask '{0}', so reverting to default".format(self.maskname))
			self.createMask(empty=True)

		# loop through all the bins, and apply the mask to the individual lightcurves
		nBins = len(self.bins)
		for i in range(nBins):
			self.bins[i].bad = self.mask[i,:]



	def fit(self, planet, star, instrument, plot=True):
		'''Take an input planet, star, and instrument; and do a fit across all wavelength bins.'''
		self.rp_over_rs = np.zeros(len(self.bins))
		self.uncertainty = np.zeros_like(self.rp_over_rs)
		self.wavelengths = np.zeros_like(self.rp_over_rs)


		#for i in range(len(self.bins)):
		#	b = self.bins[i]
		#	b.fit(planet, star, instrument, plot=plot)
		#	self.rp_over_rs[i] = b.tm.planet.rp_over_rs.value
		##	self.wavelengths[i] = b.wavelength

	def load(self, method='mcmc'):
		'''Load light curves and fits for this transmission spectrum.'''



		self.wavelengths = np.zeros(len(self.bins))
		self.fitted, self.uncertainty = {}, {}

		for i in np.arange(len(self.bins)):
			# select an individual bin
			bin = self.bins[i]
			# load its TLC and TM
			bin.load()
			# load (or create) the fit
			if method=='mcmc':
				bin.tm.slowfit(plot=False, remake=False)
			else:
				bin.tm.fastfit(plot=False, remake=False)

			# loop over the parameters, and store them in the spectrum
			for parameter in bin.tm.fitter.pdf.parameters:
				key = parameter.name
				try:
					self.fitted[key][i] = parameter.value
					self.uncertainty[key][i] = parameter.uncertainty
				except KeyError:
					self.fitted[key] = np.full(len(self.bins), None)
					self.uncertainty[key] = np.full(len(self.bins), None)
					self.fitted[key][i] = parameter.value
					self.uncertainty[key][i] = parameter.uncertainty

			self.wavelengths[i] = bin.wavelength



	def setupInitialConditions(self):

		self.initial = {}
		self.initial['planet'] = Planet(J=0.0, \
						rp_over_rs=0.107509268533, \
						rs_over_a =0.136854018274, \
						b =0.143228040337, \
						q=0.0, \
						period=3.95023867775, \
						t0=2456416.39659, \
						dt = -0.000109927092499, \
						esinw=0.0, \
						ecosw=0.0)
		self.initial['star'] = Star(u1 = 0.47, u2=0.33, temperature=6170.0, logg=4.27, metallicity=0.26)
		self.initial['instrument'] = Instrument(self.bins[0].tlc, order=2)

		# initialize a limb-darkening object
		self.ld = LD.LD(temperature=self.initial['star'].temperature.value, gravity=self.initial['star'].logg.value, metallicity=self.initial['star'].metallicity.value, directory = self.binningdirectory, plot=True)


	def setupFit(self, label='fixedGeometry', maskname='defaultMask', remake=True):

		# set the label to the kind of fit
		self.label = label

		self.applyMask(maskname)

		# make sure some initial conditions are set
		self.setupInitialConditions()

		# pull out the initial planet, star, and instrument
		p, s, i = self.initial['planet'], self.initial['star'], self.initial['instrument']

		# modify these according to what kind of a fit we want to use
		if label == 'fixedGeometry':

			# float the radius ratio
			p.rp_over_rs.float(limits=[0.05, 0.15])

			# a constant baseline
			i.C.float(value=1.0,limits=[0.9, 1.1])

			# instrument rotator angle (seems to matter!)
			i.rotatore_tothe1.float(value=0.002, limits=[-0.005, 0.005])

			# width of the whole spectrum, and the width in the wavelength range
			i.width_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])
			i.dwidth_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])

			# cross-dispersion centroid of the whole spectrum, and in the wavelength range
			i.centroid_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])
			i.dcentroid_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])

			# applied wavelength offset
			i.shift_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])

			# the sky brightness in
			i.sky_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])
			i.peak_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])


			#i.width_comparison01_tothe1.float(value=0.002, limits=[-0.005, 0.005])
			#i.sky_comparison01_tothe1.float(value=0.002, limits=[-0.005, 0.005])
			#i.centroid_comparison01_tothe1.float(value=0.002, limits=[-0.005, 0.005])
			#i.width_comparison01_tothe2.float(value=0.002, limits=[-0.005, 0.005])
			#i.width_comparison01_tothe1.float(value=0.002, limits=[-0.005, 0.005])

			# allow the limbdarkening to float [a prior for each bin will be set later]
			s.u1.float(value=s.u1.value, limits=[0.0, 1.0])
			s.u2.float(value=s.u2.value, limits=[0.0, 1.0])
			return

		if label == 'floatingGeometry':
			p.rs_over_a.float(value=0.14, limits=[0.0,1.0], shrink=1000.0)
			p.rp_over_rs.float(limits=[0.05, 0.15])
			p.b.float(value=0.8, limits=[0.0, 1.0], shrink=1000.0)
			p.dt.float(limits=np.array([-0.01, 0.01]))
			i.C.float(value=1.0,limits=[0.9, 1.1])
			i.airmass_tothe1.float(value=0.002, limits=[-0.005, 0.005])
			i.rotatore_tothe1.float(value=0.002, limits=[-0.005, 0.005])
			i.width_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])
			i.centroid_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])
			i.shift_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])
			i.dwidth_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])
			i.dcentroid_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])
			i.sky_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])
			i.peak_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])
			return

		assert(False)


	def psi(self):
		"""Quick wrapper to spit out a tuple of the initial conditions."""
		try:
			return self.initial['planet'], self.initial['star'], self.initial['instrument']
		except AttributeError:
			self.setupFit()
			return self.initial['planet'], self.initial['star'], self.initial['instrument']

	def fitBins(self, label='fixedGeometry', maskname='defaultMask', remake=False, slow=False, plot=False, **kw):
		self.speak('about to fit {0} bins with:')
		for k in locals().keys():
			self.speak('   {0} = {1}'.format(k, locals()[k]))
		self.input('are you okay with that?')
		assert(self.label == label)
		for b in self.bins:
			b.fit(plot=plot, slow=slow, remake=remake, label=label, maskname=maskname, **kw)
			assert(self.label == label)

def fastfit(inputs):
	i, kw = inputs
	t = TransmissionSpectrum(**kw)
	t.speak('starting fit for bin {0}'.format(i))
	t.bins[i].fit(plot=False, slow=False, interactive=False, **kw)
	return 'done!'

def slowfit(inputs):
	i, kw = inputs
	t = TransmissionSpectrum(**kw)
	t.speak('starting fit for bin {0}'.format(i))
	t.bins[i].fit( plot=False, slow=True, interactive=False, nburnin=500, ninference=500, **kw)
	return "done!"
