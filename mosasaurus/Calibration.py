from .imports import *
from .CCD import CCD

class Calibration(Talker):
	'''
	Calibrations are objects that store calibration data,
	including both afternoon exposures (biases, darks, flats)
	and some on-sky exposures (direct images, master science images).

	This can be thought of as bookmarked pages of
	the reducing mosasaurus' reference books,
	keeping track of calibrations that might be useful.
	'''

	def __init__(self, reducer, **kwargs):
		'''Initialize calibration object.'''

		# decide whether or not this Reducer is chatty
		Talker.__init__(self, **kwargs)

		self.speak('setting up calibrator')

		# connect to useful componenets
		self.reducer = reducer
		self.obs = self.reducer.obs
		self.display = self.reducer.display

		# create a CCD object associated with this calibration
		self.ccd = CCD(self.obs, calib=self)

		# should we be visualizing the steps?
		self.visualize = self.reducer.visualize

		# create a dictionary of full-frame calibration-relevant images
		self.images = {}

		# make all the data we need
		self.setup()

	def setup(self):
		'''Create ingredients we'll need for calibrating all images.'''

		# create master images
		self.createMasterImages()

		# figure out which are the bad pixels
		# self.createBadPixelMask()

		self.speak('calibration data are processed and ready for use')

	def createMasterImages(self, remake=False):
		'''
		Combine individual exposures into master frames
		for the various calibrations and references.
		'''

		# pull out the list of image types for which we want masters
		keys = self.obs.fileprefixes.keys()
		self.speak('creating master images for {}'.format(keys))

		# loop over all possible image types
		for k in keys:
			# create a stacked master image for everything except the science images
			#if 'science' not in k:
			self.createMasterImage(k, remake=remake)


	def createMasterImage(self, imageType=None, remake=False):
		'''
		Create a master image from a stack of them,
		using several different methods, depending on the
		type of image being processed.
		'''

		self.speak('creating a stacked master image for {}'.format(imageType))

		# if no name included, do nothing
		if imageType is None:
			self.speak('  (no image type defined; doing nothing)')
			return

		# set the CCD to a particular image type
		self.ccd.set(exposureprefix=None, imageType=imageType)

		# make sure this is a master image that is possible to make
		assert(imageType in self.obs.fileprefixes.keys())


		# we're going to save a master stacked image, and its standard deviation
		self.speak("populating the master {0} image".format(imageType))
		masterFilePrefix = os.path.join(self.obs.directory, "master_{0}".format(imageType))
		noisestring = 'StdDev'
		try:
			# has the master image already been created?
			self.images[imageType] = readFitsData(masterFilePrefix + '.fits')
			# has the master standard deviation image already been created?
			self.images[imageType+noisestring] = readFitsData(masterFilePrefix + noisestring + '.fits')
			# if both are true, then we're all set
			self.speak( "loaded {0} from {1}.fits".format(imageType, masterFilePrefix))
		except IOError:
			# create a stacked image from the appropriate file prefixes
			self.speak("creating from images " + truncate(str(self.obs.fileprefixes[imageType]), n=30))
			self.images[imageType], self.images[imageType+noisestring] = self.createStackedImage(self.obs.fileprefixes[imageType], imageType=imageType)

			# write these out to FITS files, so we can look at them in ds9
			writeFitsData(self.images[imageType], masterFilePrefix  + '.fits')
			writeFitsData(self.images[imageType+noisestring],masterFilePrefix + noisestring + '.fits')

			### FIX ME ### -- make sure the displays work nicely, for making images and movies
			#self.display.one(self.images[imageType+noisestring], clobber=True)
			#self.display.one(self.images[imageType], clobber=False)
			#self.display.single()
			#self.display.zoom()
			#self.display.scale('log', limits=[0,np.percentile(self.images[imageType],99)])
			#assert('n' not in self.input("Do you like master image {0}? [Y,n]".format(imageType)).lower())

	def createStackedImage(self, n, imageType=None, threshold=5.0, truncation=100):
		'''
		Take an outlier-rejected stack of a series of images
		(requires enough memory to hold them all).

		n is an array of fileprefixes
		imageType is a string describing the image type
		threshold is how many sigma for clipping in the stacked image
		truncation limits the number of images to be included (for large cubes)
		'''

		# if there are more than "truncation" images, take only some fraction of them
		stride = np.int(np.maximum(len(n)/truncation, 1))
		if stride > 1:
			self.speak('stacking {0}/{2} {1} images'.format(len(n),imageType,truncation))
		else:
			self.speak('stacking {0} {1} images'.format(len(n), imageType))

		# create a 3D array of images
		array = self.ccd.loadImages(n[::stride], imageType=imageType)

		# if there's only one image, simply return that image (with no noise)
		if len(array.shape) <=2:
			return array, array*0

		# calculate the outlier-rejected mean, and the 1.48*MAD for the cube
		mean, noise = zachopy.twod.stack(array, axis=0, threshold=threshold)

		#self.ccd.display.many(array, depth=0, clobber=True)
		return mean, noise

	def createBadPixelMask(self, visualize=True):
		'''Try to estimate bad pixels from a flat image. KLUDGE'''

		### STILL NEEDS COMMENTING! ###

		self.speak("populating bad pixel mask")
		badPixelFilename = self.obs.instrument.workingDirectory + 'master_BadPixels.fits'
		try:
			self.images['BadPixels'] = readFitsData(badPixelFilename)
			self.speak( "loaded bad pixel mask from {0}".format(badPixelFilename))
		except:
			self.speak( "creating bad pixel mask from the master flat frames")
			c = self.ccd#CCD(self.obs, calib=self)

			cube = []
			for n in self.obs.nWideFlat:
				c.set(n, 'WideFlat')
				cube.append(c.readData())

			cube = np.array(cube)

			median = np.median(cube,0)
			noise = np.median(np.abs(cube - median.reshape(1,cube.shape[1], cube.shape[2])), 0)
			plt.figure('bad pixel mask')
			ax = plt.subplot()
			ax.plot(median.flatten(), noise.flatten(), color='black', alpha=0.5, marker='o', markersize=4, markeredgewidth=0, linewidth=0)
			ax.set_yscale('log')
			bad = (noise < 0.05*np.sqrt(median)) | (noise == 0) | (median == 0) | (median < 0) | (self.bias() > 10000) | (self.dark() > 100)
			ax.plot(median[bad].flatten(), noise[bad].flatten(), color='red', alpha=0.5, marker='o', markersize=10, markeredgecolor='red', linewidth=0)
			ax.set_xlabel('Fluence')
			ax.set_ylabel('RMS')
			self.images['BadPixels'] = bad.astype(np.int)
			if visualize:
				self.display.one(self.images['BadPixels'])
				answer = self.input("Does the bad pixel mask seem reasonable? [Y,n]").lower()
			assert('n' not in answer)
			writeFitsData(self.images['BadPixels'], badPixelFilename)

	def bias(self):
		try:
			return self.images['bias']
		except KeyError:
			self.createMasterImage('bias')
			return self.images['bias']

	def dark(self):
		try:
			return self.images['dark']
		except KeyError:
			self.createMasterImage('dark')
			return self.images['dark']

	def science(self):
		'''
		Return the stacked science image (remaking it if necessary).
		'''
		try:
		  return self.images['science']
		except KeyError:
		  self.createMasterImage('science')
		  return self.images['science']

	def wideflat(self):
		'''
		Return the spectroscopic flat (remaking it if necessary).
		'''
		try:
			return self.images['WideFlat']
		except KeyError:
			self.createMasterImage('WideFlat')
			return self.images['WideFlat']
