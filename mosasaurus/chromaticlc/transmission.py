from Observation import Observation
import Cube
from  transit import *
import limbdarkening.phoenix as LD
from .imports import *


class TS(Talker):
	""" A transmission spectrum object, which can contain depths + uncertainties, lightcurves, covariance matrices."""

	def __init__(self, obs=None, binsize=100, remake=False, noiseassumedforplotting = 0.001, label=''):
		'''Initialize a transmission spectrum object, using a normal Observation file.'''
		Talker.__init__(self, line=200)
		self.remake= remake

		# load the observation structure for this file
		if type(obs) == str:
			self.obs = Observation(obs, nods9=True)
		else:
			self.obs = obs

		# load the
		self.name = self.obs.name
		self.binsize = binsize

		self.label = 'defaultFit'
		self.maskname = 'defaultMask'


		self.unit = 10
		self.unitstring = 'nm'
		self.noiseassumedforplotting = noiseassumedforplotting

		self.loadLCs()

	@property
	def binningdirectory(self):
		bd =  self.obs.extractionDirectory + "binby{binsize:05.0f}/".format(binsize=self.binsize)
		zachopy.utils.mkdir(bd)
		return bd

	@property
	def fitdirectory(self):
		fd = self.binningdirectory + '{label}/'.format(label=self.label)
		zachopy.utils.mkdir(fd)
		return fd

	@property
	def maskdirectory(self):
		md = self.fitdirectory + "{maskname}/".format(maskname=self.maskname)
		zachopy.utils.mkdir(md)
		return md

	def loadLCs(self):
		try:
			# first try loading light curves that have been saved into the directory structure
			possibleTLCs = glob.glob(self.binningdirectory + '*to*/TLC.npy')
			assert(len(possibleTLCs) > 0)
			chunks = []
			for file in possibleTLCs:
				chunks.append(file.split('/')[-2])

		except:
			lcDirectory = self.obs.extractionDirectory + 'lc_binby{0}'.format(self.binsize) + '/'
			possibleTLCs = glob.glob(lcDirectory + '*.lightcurve')
			assert(len(possibleTLCs) > 0)
			chunks = []
			for file in possibleTLCs:
				chunks.append((file.split('/')[-1].split('.lightcurve')[0]).split('_')[-1])

		bins = []
		wavelengths = []
		for trimmed in chunks:
			left = np.int(trimmed.split('to')[0])
			right = np.int(trimmed.split('to')[-1])
			self.speak("creating a bin between {left} and {right}".format(left=left, right=right))
			try:
				bins.append(Bin(self, left=left, right=right))
				wavelengths.append(bins[-1].wavelength)
			except IOError:
				pass
		self.bins = np.array(bins)[np.argsort(wavelengths)]


	def toArray(self, key):
		'''Create an image array of a TLC variable.'''
		nBins = len(self.bins)
		nTimes = len(self.bins[0].tlc.flux)
		a = np.zeros((nBins, nTimes))
		for i in range(nBins):
			tlc = self.bins[i].tlc
			try:
				a[i,:] = tlc.__dict__[key]
			except:
				a[i,:] = tlc.externalvariables[key]
		return a

	def createMask(self):

		a = raw_input("What would you like to call this mask?\n ")
		maskname = a.replace(" ","")

		keys = ['peak_', 'sky_', 'width_', 'centroid_']
		# one column for light curves, one for target, one for comparison
		nColumns = 3
		nRows = len(keys)
		plt.figure('masking')
		ip = zachopy.displays.iplot.iplot(nRows, nColumns)

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
				self.speak(" in bin {0}, {1} points exceeded {2} x {3}".format(i, np.sum(bad), threshold, noise))
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
					self.speak("  Click at the two corners of a box you'd like to mask.\n")
					clicks = ip.getMouseClicks(2)

					rows = (clicks[0].ydata, clicks[1].ydata)
					top = np.int(np.max(rows))
					bottom = np.maximum(np.int(np.min(rows)), 0)

					columns =  (clicks[0].xdata, clicks[1].xdata)
					right = np.int(np.max(columns))
					left = np.maximum(np.int(np.min(columns)), 0)


					mask[bottom:top, left:right] = mask[bottom:top, left:right] | self.bins[0].tlc.flags['custom']

				elif 's' in answer:
					self.speak("  Click at the two corners of a box you'd like to unmask.\n")
					clicks = ip.getMouseClicks(2)

					rows = np.round((clicks[0].ydata, clicks[1].ydata))
					top = np.int(np.max(rows))
					bottom = np.maximum(np.int(np.min(rows)), 0)

					columns =  np.round((clicks[0].xdata, clicks[1].xdata))
					right = np.int(np.max(columns))
					left = np.maximum(np.int(np.min(columns)), 0)

					mask[bottom:top, left:right] -= mask[bottom:top, left:right] & self.bins[0].tlc.flags['custom']

				elif 'r' in answer:
					self.speak("Okay, refitting. It may take a while!")
					self.mask = mask
					self.applyMask(maskname=maskname)
					self.fitRigid()
				elif 'f' in answer:
					keepgoing = False
				else:
					unreasonableanswer = True
					self.speak("  I'm sorry, I didn't quite understand that.")

		self.mask = mask
		self.applyMask(maskname=maskname)

	def applyMask(self, maskname=None):
		if maskname is not None:
			self.maskname = maskname
		np.save(self.maskdirectory + 'mask.npy', self.mask)
		nBins = len(self.bins)
		for i in range(nBins):
			self.bins[i].tlc.bad = self.mask[i,:]

	def fit(self, planet, star, instrument, plot=True):
		'''Take an input planet, star, and instrument; and do a fit across all wavelength bins.'''
		self.rp_over_rs = np.zeros(len(self.bins))
		self.uncertainty = np.zeros_like(self.rp_over_rs)
		self.wavelengths = np.zeros_like(self.rp_over_rs)
		self.ld = LD.LD(temperature=star.temperature.value, gravity=star.logg.value, metallicity=star.metallicity.value, directory = self.binningdirectory, plot=True)

		for i in range(len(self.bins)):
			b = self.bins[i]
			b.fit(planet, star, instrument, plot=plot)
			self.rp_over_rs[i] = b.tm.planet.rp_over_rs.value
			self.uncertainty[i] =  b.tm.planet.rp_over_rs.uncertainty
			self.wavelengths[i] = b.wavelength


			plt.savefig(b.fittingdirectory + 'lightcurve.pdf')
			self.speak(b.fittingdirectory + 'lightcurve.pdf')


	def load(self):
		'''Load light curves and fits for this transmission spectrum.'''
		self.rp_over_rs = np.zeros(len(self.bins))
		self.uncertainty = np.zeros_like(self.rp_over_rs)
		self.wavelengths = np.zeros_like(self.rp_over_rs)

		for i in np.arange(len(self.bins)):
			# select this bin
			bin = self.bins[i]
			bin.load()
			self.rp_over_rs[i] = bin.tm.planet.rp_over_rs.value
			self.uncertainty[i] =  bin.tm.planet.rp_over_rs.uncertainty
			self.wavelengths[i] = bin.wavelength


	def setupSuperPlot(self):
		# set up grid of plots
		gs = plt.matplotlib.gridspec.GridSpec(3, 1,wspace=0,hspace=0.05,height_ratios=[1,2,2])

		# create the plot for the spectrum of the star
		self.ax_spectrum = plt.subplot(gs[0])
		self.ax_spectrum.set_ylabel('Raw Detected Flux (photons/nm/exposure)')
		plt.setp(self.ax_spectrum.get_xticklabels(), visible=False)

		# create the plot for the light curves
		self.ax_lc = plt.subplot(gs[1], sharex=self.ax_spectrum)
		self.ax_lc.set_ylabel('Time from Mid-Transit\n(hours)')
		plt.setp(self.ax_lc.get_xticklabels(), visible=False)

		# create the plot for the transmission spectrum
		self.ax_ts = plt.subplot(gs[2], sharex=self.ax_spectrum)
		self.ax_ts.set_xlim((self.bins[0].left - self.binsize/2)/self.unit, (self.bins[-1].right +  self.binsize/2)/self.unit)
		self.ax_ts.set_xlabel('Wavelength ({0})'.format(self.unitstring))
		self.ax_ts.set_ylabel('Transit Depth (%)')



	def plot(self):
		self.setupSuperPlot()
		colors = []
		def normalize(flux):
			ratio = 0.8
			one = (flux-1.0)/np.mean(self.rp_over_rs)**2
			return self.binsize/self.unit*ratio*(one+0.5)

		# plot the spectrum
		try:
			wavelength, spectrum = np.load(self.obs.extractionDirectory + 'medianSpectrum.npy')
			self.ax_spectrum.plot(wavelength/self.unit, spectrum*self.unit, linewidth=3, alpha=0.5, color='black')
		except:
			self.speak('no median spectrum could be found')
		for i in np.arange(len(self.bins)):
			# select this bin
			bin = self.bins[i]

			# plot the model for this bin
			#time, planetmodel, instrumentmodel = bin.tlc.TM.smooth_model()
			#kw = {'color':zachopy.color.nm2rgb([bin.left/self.unit, bin.right/self.unit], intensity=0.25), 'linewidth':3, 'alpha':0.5}
			#self.ax_lc.plot(normalize(planetmodel) + bin.wavelength/self.unit, bin.tlc.TM.planet.timefrommidtransit(time), **kw)

			# plot the (instrument corrected) datapoints
			kw = {'marker':'.', 'color':zachopy.color.nm2rgb([bin.left/self.unit, bin.right/self.unit]), 'alpha':0.25, 'linewidth':0, 'marker':'o'}
			self.ax_lc.plot(normalize(bin.tlc.corrected()[bin.tlc.bad == False]) + bin.wavelength/self.unit, bin.tlc.timefrommidtransit()[bin.tlc.bad == False], **kw)
			colors.append(kw['color'])

			self.speak(bin)
			self.speak(bin.tlc.TM.planet)

			width = 3
			self.ax_ts.errorbar(self.wavelengths[i]/self.unit, self.rp_over_rs[i], self.uncertainty[i], marker='o', color=zachopy.color.nm2rgb([bin.left/self.unit, bin.right/self.unit]), markersize=10, linewidth=width, elinewidth=width, capsize=5, capthick=width)

	def fitRigid(self, plot=False):
		plt.ion()
		self.label = 'fixedGeometry'

		p = Planet(J=0.0, \
						rp_over_rs=0.107509268533, \
						rs_over_a =0.136854018274, \
						b =0.143228040337, \
						q=0.0, \
						period=3.95023867775, \
						t0=2456416.39659, \
						dt = -0.000109927092499, \
						esinw=0.0, \
						ecosw=0.0)
		s = Star(u1 = 0.47, u2=0.33, temperature=6170.0, logg=4.27, metallicity=0.26)
		i = Instrument(self.bins[0].tlc, order=2)

		p.rp_over_rs.float(limits=[0.05, 0.15])

		# a constant baseline
		i.C.float(value=1.0,limits=[0.9, 1.1])
		#i.airmass_tothe1.float(value=0.002, limits=[-0.005, 0.005])

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

		s.u1.float(value=s.u1.value, limits=[0.0, 1.0])
		s.u2.float(value=s.u2.value, limits=[0.0, 1.0])


		self.fit(p, s, i, plot=plot)

	def fitFlexible(self, filename='wasp94_140805.obs',binsize=100, remake=False, plot=False):
		plt.ion()
		self.label = 'floatingGeometry'
		p = Planet(J=0.0, \
						rp_over_rs=0.10745, \
						rs_over_a=0.136146284482, \
						b = 0.171287057989, \
						q=0.0, \
						period=3.95023867775, \
						t0=2456416.39659, \
						esinw=0.0, \
						ecosw=0.0)
		s = Star(u1 = 0.47, u2=0.33, temperature=6170.0, logg=4.27, metallicity=0.26)
		i = Instrument(self.bins[0].tlc, order=2)

		p.rs_over_a.float(limits=[0.0,1.0])
		p.rp_over_rs.float(limits=[0.05, 0.15])
		p.b.float(limits=[0.0, 1.0])
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
		#i.width_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])
		#i.width_target_tothe2.float(value=0.002, limits=[-0.005, 0.005])
		#i.centroid_tothe1.float(value=0.002, limits=[-0.005, 0.005])
		#i.centroid_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])
		#i.centroid_target_tothe2.float(value=0.002, limits=[-0.005, 0.005])
		#i.sky_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])
		#i.sky_target_tothe2.float(value=0.002, limits=[-0.005, 0.005])
		#i.peak_target_tothe1.float(value=0.002, limits=[-0.005, 0.005])
		#i.width_comparison01_tothe1.float(value=0.002, limits=[-0.005, 0.005])
		#i.sky_comparison01_tothe1.float(value=0.002, limits=[-0.005, 0.005])
		#i.centroid_comparison01_tothe1.float(value=0.002, limits=[-0.005, 0.005])
		#i.width_comparison01_tothe2.float(value=0.002, limits=[-0.005, 0.005])
		#i.width_comparison01_tothe1.float(value=0.002, limits=[-0.005, 0.005])

		self.fit(p, s, i, plot=plot)



def load(filename='wasp94_140805.obs',binsize=500, label='fixedGeometry'):
	ts = TS(filename, binsize=binsize, label=label)
	ts.label = label
	ts.loadLCs()
	#ts.load()
	return ts

def plotTS(filename='wasp94_140805.obs',binsize=500):
	ts = TS(filename, binsize=binsize, label='fixedGeometry')
	ts.loadLCs()
	ts.load()
	ts.plot()

def determineParameters(filename='wasp94_140805.obs',binsize=500):
	ts = TS(filename, binsize=binsize, label='floatingGeometry')
	ts.loadLCs()
	ts.label = 'floatingGeometry'
	ts.load()
	ts.plot()

	keys = ['rs_over_a', 'b', 'dt']
	n = len(keys)
	dict = {}
	for k in keys:
		dict[k+"_value"] = [bin.tm.planet.__dict__[k].value for bin in ts.bins]
		dict[k+"_uncertainty"] = [bin.tm.planet.__dict__[k].uncertainty for bin in ts.bins]

	self.speak("")
	self.speak("The median values of the {binsize} angstrom fits for {filename} are:".format(binsize=binsize, filename=filename))


	plt.figure('geometric parameters')
	gs = plt.matplotlib.gridspec.GridSpec(n,n)
	for i in range(n):
		for j in range(n):
			ax = plt.subplot(gs[j,i])
			x = dict[keys[i] + '_value']
			y = dict[keys[j] + '_value']
			xerr = dict[keys[i] + '_uncertainty']
			yerr = dict[keys[j] + '_uncertainty']
			xstd = np.std(x)
			ystd = np.std(y)
			xmed = np.median(x)
			ymed = np.median(y)
			nsigma=3
			ax.errorbar(x, y, xerr=xerr, yerr=yerr, marker='o', linewidth=0, elinewidth=3, capthick=3, color='black', alpha=0.3)
			ax.plot(xmed, ymed, marker='o', markersize=20, alpha=0.5, color='red')
			ax.set_xlabel(keys[i])
			ax.set_ylabel(keys[j])
			ax.set_xlim(xmed - nsigma*xstd, xmed + nsigma*xstd)
			ax.set_ylim(ymed - nsigma*ystd, ymed + nsigma*ystd)
		self.speak("{0:>20} = {1}".format(keys[i], xmed))
