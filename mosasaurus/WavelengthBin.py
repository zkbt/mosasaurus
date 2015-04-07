

class WavelengthBin(Talker):
	'''Bin object to store everything for one transmission spectrum bin.'''
	def __init__(self, ts, left=None, right=None):
		'''Initialize a Bin object.'''

		Talker.__init__(self)

		# a Bin has wavelength definitions
		self.left = left
		self.right = right
		assert(self.left is not None)

		# a Bin is associated with a parent
		self.TS = ts
		self.unit = self.TS.unit
		self.unitstring = self.TS.unitstring

		# give an update
		self.speak("Initializing a spectroscopic bin covering {left} to {right} {unitstring}, storing data in {datadirectory} and fits in {fittingdirectory}".format(left=self.left/self.unit, right=self.right/self.unit, unitstring=self.unitstring, datadirectory=self.datadirectory, fittingdirectory=self.fittingdirectory))

		# load the light curves, either from processed file or from a raw file somewhere
		try:
			# try to load the processed light curve file
			assert(self.TS.remake == False)
			self.readProcessedLC()
		except:
			self.readRawLC()

		# load a model, either a precalculated one, or a newly generated one
		try:
			self.speak('attempting to load a transit model from {0}'.format(self.fittingdirectory))
			self.tm = TM(directory=self.fittingdirectory)
			self.tlc.linkModel(self.tm)
		except:
			self.speak("could not load a transit model from {0}".format( self.fittingdirectory))
			pass

	@property
	def wavelength(self):
		return (self.right + self.left)/2.0

	@property
	def binsize(self):
		return self.right - self.left


	def readProcessedLC(self):
		self.speak( "attempting to load a processed light curve from {0}".format( self.datadirectory))
		self.tlc = TLC(directory=self.datadirectory, name=self.TS.name, left=self.left, right=self.right)
		assert(self.tlc is not None)
		self.speak( "   ...success!")

	def readRawLC(self):#, filtersize=3):
		lcFile = self.TS.obs.extractionDirectory + 'lc_binby{0}'.format(self.binsize) + '/' + "{0:05.0f}to{1:05.0f}.lightcurve".format(self.left, self.right)
		self.speak("attempting to load a raw light curve from {0}".format(lcFile))
		table = astropy.io.ascii.read(lcFile)
		arrays = {}
		for k in table.colnames:
			if 'col' not in k:
				arrays[k] = table[k].filled().data
		self.speak('before')
		self.tlc = TLC(name=self.TS.name, left=self.left, right=self.right, **arrays)
		self.speak('after')
		self.speak('...success!')
		self.tlc.save(self.datadirectory)

	@property
	def datadirectory(self):
		dd = self.TS.binningdirectory + "{0:05.0f}to{1:05.0f}/".format(self.left, self.right)
		zachopy.utils.mkdir(dd)
		return dd

	@property
	def fittingdirectory(self):
		md = self.TS.maskdirectory + "{0:05.0f}to{1:05.0f}/".format(self.left, self.right)
		zachopy.utils.mkdir(md)
		return md


	def fit(self, planet, star, instrument, plot=True):
		print "  Fitting light curve covering {left} to {right} {unitstring}.".format(left=self.left/self.unit, right=self.right/self.unit, unitstring=self.unitstring)
		plt.ion()
		self.planet = copy.deepcopy(planet)
		self.depthassumedforplotting = self.planet.rp_over_rs.value**2
		self.star = copy.deepcopy(star)

		u1, u2, du1dt, du2dt = self.TS.ld.quadratic(self.left, self.right, nudge=True)
		self.star.u1.value = u1
		self.star.u2.value = u2
		self.star.u1.uncertainty = np.abs(du1dt)*200.0
		self.star.u2.uncertainty = np.abs(du2dt)*200.0

		self.instrument = copy.deepcopy(instrument)
		self.tm = TM(self.planet, self.star, self.instrument, depthassumedforplotting=self.depthassumedforplotting , directory=self.fittingdirectory)
		self.tlc.linkModel(self.tm)
		self.tm.fastfit(plot=plot)
		#if plot:
		#ls	self.tlc.DiagnosticsPlots()
		self.save()

	def save(self):
		self.tlc.save(self.datadirectory)
		self.tm.save(self.fittingdirectory)

	def load(self):
		self.tlc = TLC(directory=self.datadirectory)
		self.tm = TM(directory=self.fittingdirectory)
		self.tlc.linkModel(self.tm)

	def __str__(self):
		s = "spectroscopic bin:\n {left:.0f}-{right:.0f}nm\n {n} light curve points\n".format(left=self.left/self.unit, right = self.right/self.unit, n=len(self.tlc.bjd))
		try:
			for key in self.tm.planet.__dict__.keys():
				p = self.tm.planet.__dict__[key]
				if p.uncertainty > 0:
					s += "     {key} = {value}\pm{uncertainty}\n".format(key=key, value=p.value, uncertainty=p.uncertainty)
		except:
			"     no model defined."
		return s
