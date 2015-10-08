from imports import *
from TransmissionSpectrum import WithTLCs, TransmissionSpectrum
from WavelengthBin import WavelengthBin
import transit

class CombinedObs(Talker):

	def __init__(self, listofobs):
		obs = listofobs[0]
		self.name = obs.name
		self.night = 'combined'
		self.baseDirectory = obs.baseDirectory
		self.workingDirectory = obs.baseDirectory + 'working/combinationof{0}observations/'.format(len(listofobs))
		zachopy.utils.mkdir(self.workingDirectory)
		self.extractionDirectory = self.workingDirectory + obs.extractionDirectory.split('/')[-2] + '/'
		zachopy.utils.mkdir(self.extractionDirectory)

class CombinedTransmissionSpectrum(TransmissionSpectrum):
	def __init__(self, method='lm', label='fixedGeometry', maskname='defaultMask', files=['wasp94_140801.obs', 'wasp94_140805.obs', 'wasp94_140809.obs'], binsize=100):
		Talker.__init__(self)


		# keep track of the files
		self.files = files

		# populate list of transmission spectra
		self.spectra = []


		# set the mask and label names
		self.maskname = maskname
		self.label = label

		# loop over spectra, loading each
		self.speak('creating a combined transmission spectrum from files:')
		for f in files:
			self.speak('  {0}'.format(f))
			spectrum = WithTLCs(f, binsize=binsize)
			self.spectra.append(spectrum)

		# create a night-agnostic observation, that ties them all together
		self.obs = CombinedObs([s.obs for s in self.spectra])

		# set the binsize
		self.setBinsize(self.spectra[0].binsize)


		# create wavelength bins to cover all options among constiuent spectra
		binoptions = []
		for s in self.spectra:
			for b in s.bins:
				binoptions.append(b.__repr__().split('|')[-1].split('>')[0])
				print binoptions
		bins, wavelengths = [], []
		# create the bins, but don't populate with light curves
		for b in np.unique(binoptions):
			left, right = np.array(b.split('nm')[0].split('to')).astype(np.int)*10
			bins.append(WavelengthBin(self, left=left, right=right))
			wavelengths.append(bins[-1].wavelength)
		s = np.argsort(wavelengths)
		self.bins = np.array(bins)[s]
		self.wavelengths = np.array(wavelengths)[s]

		# create a dictionary to store all the TLCs for each wavelength
		self.archiveoftlcs = {}
		for s in self.spectra:
			for b in s.bins:
				w = b.wavelength
				try:
					self.archiveoftlcs[w].append(b.tlc)
				except KeyError:
					self.archiveoftlcs[w] = [b.tlc]

	def fitMonochromatic(self, w, initialConditions=None, wobbly=True,  remake=False):
		''' fit a wavelength bin, across all observations, given some initial conditions '''

		# create an empty list of tlcs
		tlcs = []

		if wobbly:
			self.label = 'floatingGeometry'
		else:
			self.label = 'fixedGeometry'

		# loop over the tlcs, and create a model for each
		for orig in self.archiveoftlcs[w]:

			# define the input objects
			planet = transit.Planet(**initialConditions.planetkw)

			synthesizerdirectory = '{0}{1}/'.format(self.fitdirectory,self.bins[self.w2bin(w)][0].identifier)
			zachopy.utils.mkdir(synthesizerdirectory)
			
			# assign an epoch to the TLC
			tlc = orig.splitIntoEpochs(planet,
				newdirectory=synthesizerdirectory)[0]
			tlcs.append(tlc)

			# float the planetary parameters
			planet.k.float(limits=[0.05, 0.15])
			planet.b.float(limits=[0.0, 1.0])
			planet.rsovera.float(limits=[0.01, 0.5])
			planet.period.float(limits=[planet.period.value-1e-5, planet.period.value+1e-5])

			instrument = transit.Instrument(tlc=tlc, **initialConditions.instrumentkw)

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
			instrument.shift_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])


			# the sky brightness in
			instrument.sky_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])
			instrument.peak_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])

			# pull out the limbdarkening coefficients
			u1, u2 = initialConditions.ld.quadratic(tlc.left, tlc.right)
			star = transit.Star(u1=u1, u2=u2)

			# create the transit model
			tm = transit.TM(planet, star, instrument, directory=tlc.directory)
			tm.linkLightCurve(tlc)

		self.synthesizer = transit.LM(tlcs=tlcs, directory=synthesizerdirectory)
		for thing in ['period', 't0',  'b', 'rsovera', 'semiamplitude']:
			self.synthesizer.tieAcrossEpochs(thing)
			self.synthesizer.tieAcrossTelescopes(thing)

		for thing in ['k', 'u1', 'u2']:
			self.synthesizer.tieAcrossEpochs(thing)

		self.synthesizer.speak('the starting parameters are')
		self.synthesizer.printParameters()

		try:
			self.archiveofpdfs
		except AttributeError:
			self.archiveofpdfs = {}

		self.synthesizer.fit(remake=remake)
		self.archiveofpdfs[w] = self.synthesizer.pdf


	def load(self, method='lm'):
		self.speak('loading all spectra')
		for s in self.spectra:
			s.load(method)

	def combine(self):

		self.fitted, self.uncertainty, self.chisq, self.dof = {}, {}, {}, {}

		# loop over wavelengths

		for i in range(len(self.wavelengths)):
			w = self.wavelengths[i]
			self.speak('combining all available spectra for {0}'.format(w))
			f2combine, u2combine = {},{}

			for s in self.spectra:
				ok = (s.wavelengths == w).nonzero()
				self.speak('wavelength {0} is bin {1}'.format(w, ok))
				assert(len(ok) == 1)
				for k in s.fitted.keys():
					try:
						f2combine[k], u2combine[k]
					except KeyError:
						f2combine[k], u2combine[k] = [], []
					f2combine[k].append(s.fitted[k][ok][0])
					u2combine[k].append(s.uncertainty[k][ok][0])

			for k in f2combine.keys():
				f2combine[k] = np.array(f2combine[k])
				u2combine[k] = np.array(u2combine[k])
				try:
					self.fitted[k], self.uncertainty[k]
				except:
					self.fitted[k] = np.empty(len(self.wavelengths))
					self.uncertainty[k] = np.empty(len(self.wavelengths))
					self.chisq[k] = np.empty(len(self.wavelengths))
					self.dof[k] = np.empty(len(self.wavelengths))

				self.fitted[k][i] = np.sum(f2combine[k]/u2combine[k]**2)/np.sum(1.0/u2combine[k]**2)
				self.chisq[k][i] = np.sum((f2combine[k] - self.fitted[k][i])**2/u2combine[k]**2)
				self.dof[k][i] = max(len(f2combine[k]) - 1.0, 1)
				rescaling = max(np.sqrt(self.chisq[k][i]/self.dof[k][i]), 1)

				self.uncertainty[k][i] = 1.0/np.sqrt(np.sum(1.0/u2combine[k]**2))*rescaling

			print f2combine
			print u2combine
			print

				######## PICK UP HERE!

		pass

	def table(self):
		form = '{0:>20}{1:>20}{2:>20}{3:>25}'
		print form.format('left', 'right', 'rp_over_rs', 'rp_over_rs_error')

		for i in range(len(self.bins)):
			b = self.bins[i]
			if i == 0:
				print form.format(b.unitstring, b.unitstring, 'unity', 'unity')
			print form.format(b.left/b.unit, b.right/b.unit, self.fitted['k'][i], self.uncertainty['k'][i])

	@property
	def mask(self):
		"Indicate the mask being using for all spectra."
		return self._mask

	@mask.setter
	def mask(self, value):
		"Specify the mask to use (for all spectra)."

		self.speak('setting mask for all spectra to {0}'.format(value))
		# update the maskname hidden value
		self._mask = value

		# loop over the transmission spectra
		for spectrum in self.spectra:
			self.speak('setting mask for {spectrum} to {mask}'.format(spectrum=spectrum, mask=value))
			# use the transmission spectrum's applyMask method to load the appropriate mask and apply it to the light curves
			spectrum.applyMask(value)

	@property
	def label(self):
		"Indicate the fit being used for all spectra."
		return self._label

	@label.setter
	def label(self, value):
		"Specify the mask to use (for all spectra)."

		self.speak('setting label for all spectra to {0}'.format(value))
		self._label = value
		for spectrum in self.spectra:
			self.speak('setting label for {spectrum} to {label}'.format(spectrum=spectrum, label=value))
			spectrum.label = value

	def __repr__(self):
		"""How should this object be represented (e.g. when shown as an element in a list)"""
		s = "<CombinedTransmissionSpectrum "
		for t in self.spectra:
			s += '{0}|'.format( repr(t))
		s = s.strip('|') +  '>'
		return s
