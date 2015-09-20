from imports import *
from transit import TLC, TM

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
		#self.speak("Initializing a spectroscopic bin covering {left} to {right} {unitstring}, as part of {ts}".format(left=self.left/self.unit, right=self.right/self.unit, unitstring=self.unitstring, ts=self.TS))


	def __repr__(self):
		"""How should this object be represented (e.g. when shown as an element in a list)"""

		return '<WavelengthBin {name}|{night}|{left}to{right}{unitstring}>'.format(left=self.left/self.unit, right=self.right/self.unit, unitstring=self.unitstring, name=self.TS.obs.name,night=self.TS.obs.night)

	@property
	def wavelength(self):
		return (self.right + self.left)/2.0

	@property
	def binsize(self):
		return self.right - self.left

	def readTLC(self, remake=False):
		'''read the TLC for this bin, either from the preprocessed object (fast) or a raw file (slower)'''
		try:
			# try to load the processed light curve file
			assert(self.TS.remake == False)
			assert(remake==False)
			self.readProcessedLC()
		except:
			self.readRawLC()

	def readProcessedLC(self):
		'''read a light curve that has already been loaded and saved as a TLC object'''
		self.speak( "attempting to load a processed light curve from {0}".format( self.datadirectory))
		self.tlc = TLC(directory=self.datadirectory, name=self.TS.name, left=self.left, right=self.right)
		assert(self.tlc is not None)
		self.speak( "   ...success!")

	def readRawLC(self):#, filtersize=3):
		'''read a raw light curve from a .lightcurve text file'''
		lcFile = self.TS.rawlightcurvedirectory + "{0:05.0f}to{1:05.0f}.lightcurve".format(self.left, self.right)
		self.speak("attempting to load a raw light curve from {0}".format(lcFile))
		table = astropy.io.ascii.read(lcFile)
		arrays = {}
		for k in table.colnames:
			if 'col' not in k:
				arrays[k] = table[k].data
		self.speak('before')
		self.tlc = TLC(name=self.TS.name, left=self.left, right=self.right, directory=self.datadirectory, **arrays)
		self.speak('after')
		self.speak('...success!')
		self.tlc.save(self.datadirectory)

	@property
	def datadirectory(self):
		dd = self.TS.lightcurvedirectory + "{0:05.0f}to{1:05.0f}/".format(self.left, self.right)
		zachopy.utils.mkdir(dd)
		return dd

	@property
	def fittingdirectory(self):
		md = self.TS.fitdirectory + "{0:05.0f}to{1:05.0f}/".format(self.left, self.right)
		zachopy.utils.mkdir(md)
		return md



	def fit(self, plot=True, slow=False, label='fixedGeometry', remake=True, maskname='defaultMask', **kwargs):
		"""Take an input planet, star, and instrument -- and fit the bin's light curve!"""

		# say what we're doing
		if slow:
			method='MCMC'
		else:
			method='LM'
		self.speak('fitting {0}, using {method}'.format(self, method=method))

		# mask any additional points as specified in the mask

		self.TS.setupFit(label=label, remake=remake, maskname=maskname)


		planet, star, instrument = self.TS.psi()

		#
		try:
			self.tlc
		except AttributeError:
			self.readTLC()

		# apply the spectrum-level mask
		self.tlc.bad *= self.TS.mask[self.TS.w2bin[self.wavelength], :]


		# initialize some structures
		self.planet = copy.deepcopy(planet)
		self.star = copy.deepcopy(star)
		self.instrument = copy.deepcopy(instrument)

		# set up the depth that will be assumed for making plots (so that multiple wavelength bins line up)
		self.depthassumedforplotting = self.planet.rp_over_rs.value**2


		# initialize the limb darkening values and uncertainties, with prior from atmosphere models
		u1, u2, du1dt, du2dt = self.TS.ld.quadratic(self.left, self.right, nudge=True)
		self.star.u1.value = u1
		self.star.u2.value = u2
		self.star.u1.uncertainty = np.abs(du1dt)*200.0
		self.star.u2.uncertainty = np.abs(du2dt)*200.0

		# create the transit model, using the input structures
		self.tm = TM(self.planet, self.star, self.instrument, depthassumedforplotting=self.depthassumedforplotting , directory=self.fittingdirectory)

		# link the model and the light curve
		self.tlc.linkModel(self.tm)

		if slow:
			self.tm.slowfit(plot=plot, remake=remake, **kwargs)
		else:
			self.tm.fastfit(plot=plot, remake=remake,  **kwargs)

		#if plot:
		#ls	self.tlc.DiagnosticsPlots()
		self.save()

	def save(self):
		self.tlc.save(self.datadirectory)
		self.tm.save(self.fittingdirectory)

	def readTM(self):
		self.tm = TM(directory=self.fittingdirectory)

	def load(self):
		self.readTLC()
		self.readTM()
		self.tlc.linkModel(self.tm)
