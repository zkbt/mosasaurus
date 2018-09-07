from .imports import *

class CCD(Talker):
    '''
    CCD object handles every related to an individual CCD exposure.

    You can think of this as a simulated CCD on which our dear
    mosasaurus can play back any image.
    '''

    def __init__(self, obs, exposureprefix=0, imageType=None, calib=None, **kwargs):
        '''
        Initialized the CCD object, connected to a given observation set,
        and with a name n (originally a number, now generalized to a name)
        and an imageType (which sets what basic calibrations should be applied.)'''


        # decide whether or not this Reducer is chatty
        Talker.__init__(self, **kwargs)


        # initialize the basic stuff we need
        self.obs = obs
        self.instrument = self.obs.instrument
        self.calib = calib
        self.display = self.calib.display

        # make sure the stitched directory is defined
        self.stitchedDirectory = os.path.join(self.obs.directory, 'stitched')
        mkdir(self.stitchedDirectory )


        # set up options
        self.flags = {'subtractbias':True, 'subtractdark':True, 'multiplygain':False}
        self.visualize = self.obs.reducer.visualize

        # point this CCD to right file, if provided
        self.set(exposureprefix, imageType)


    def set(self, exposureprefix=None, imageType=None):
        '''Point this CCD object at a specific file.'''

        # keep track if this is some special kind of image
        self.imageType = imageType
        self.exposureprefix = exposureprefix

        # define the file prefix
        if exposureprefix is not None:

            # define a nickname for referring to this image
            if self.imageType is None:
                self.label = 'unknown'
            else:
                self.label = self.imageType
            self.name = '{}_{}'.format(self.label, self.exposureprefix)

            # define a stitched filename
            self.stitched_filename = os.path.join(self.stitchedDirectory,'{}.fits'.format(self.name))

            # set up directories
            if self.instrument.zapcosmics:
                self.cosmicsDirectory = os.path.join(self.obs.directory, 'cosmics')
                self.cosmicsFilename =  os.path.join(self.cosmicsDirectory, 'rejectedpercolumn_{}.npy'.format(self.name))

        # empty out the header and data variables
        self.header = None
        self.data = None

    def readHeader(self, exposureprefix=None, imageType=None):
        '''Read in the header for this image.'''

        # set the CCD, if necessary
        if exposureprefix is not None:
            self.set(exposureprefix, imageType)

        # read one of the two images to get header
        filename = self.instrument.prefix2files(exposureprefix)[0]
        hdu = astropy.io.fits.open(filename)
        header = hdu[0].header
        hdu.close()

        # store the header in this object
        self.header = header

        self.speak("read header from {0}".format(self.name))

        # return the header
        return self.header

    def writeData(self):
            self.speak('saving image to {0}'.format(self.stitched_filename))
            writeFitsData(self.data, self.stitched_filename)

    def readData(self, exposureprefix=None, imageType=None):
        '''Read in the image data for this exposure (create stitched exposure if necessary.)'''

        # set the CCD, if necessary
        if exposureprefix is not None:
            self.set(exposureprefix, imageType)

        # try loading a stitched image, otherwise create one from scratch
        try:
            # does a stitched image already exist?
            self.data = readFitsData(self.stitched_filename)
        except IOError:
            # if not, you'll need to create one
            self.createStitched()

        ### FIX ME ### (come up with a better solution for cosmic ray mitigation)
        if self.instrument.zapcosmics:
            if imageType == 'science':
                self.cosmicdiagnostic = np.load(self.cosmicsFilename)

        # print status
        self.speak("read image from {0}".format(self.name))

        # return the image
        return self.data


    def rejectCosmicRays(self, remake=False, threshold=7.5, visualize=False, nBeforeAfter=5):
        '''Stitch all science images, establish a comparison noise level for each pixel.'''

        # make sure a cosmics directory exists
        cosmics_directory = os.path.join(self.obs.directory, 'cosmics/')
        mkdir(cosmics_directory)

        # figure out how many images to consider
        nImages = 2*nBeforeAfter + 1
        imageType = 'ScienceUnmitigated'
        exposureprefix = self.exposureprefix

        try:
            #print "testing"
            ok = remake == False
            #print 'remake'
            ok = ok & os.path.exists(self.stitched_filename)
            #print 'fits'
            ok = ok & os.path.exists(self.cosmicsFilename)
            #print 'rejected'
            #print ok
            assert(ok)
            self.speak('a cosmic-rejected stitched/Science{0:04.0f}.fits already exists!'.format(n))

        except (IOError,AssertionError):
            # what's the index of this exposure?
            i = self.obs.exposures['science'].loc[exposureprefix].index

            iComparison = np.arange(np.maximum(0, i-nBeforeAfter), np.minimum(i+nBeforeAfter, len(self.obs.exposures['science'])))
            prefixesComparision = self.obs.fileprefixes['science'][iComparison]

            #nComparison = np.arange(-nBeforeAfter + n, nBeforeAfter+n+1, 1)
            #nComparison = nComparison[(nComparison >= np.min(self.obs.nScience)) & (nComparison <= np.max(self.obs.nScience))]
            comparison = self.loadImages(prefixesComparision, imageType=imageType)
            shape = np.array(comparison.shape)
            axis = 0
            shape[axis] = 1

            self.speak('comparing image {0} to images {1} to remove cosmic rays'.format(exposureprefix, prefixesComparision))
            image = self.data
            #image = comparison[nComparison == n,:,:].squeeze()


            # calculate median and noise of comparisons
            med = np.median(comparison, axis=axis)
            noise = np.maximum(1.48*np.median(np.abs(comparison - med.reshape(shape)), axis=axis), 1.0)

            bad = (image - med)/noise > threshold
            corrected = image + 0.0
            corrected[bad] = med[bad]
            if visualize:
                self.display.replace(image,0)
                self.display.replace(image - corrected,1)
                self.display.replace(corrected,2)
                self.display.scale(mode='zscale')
                self.display.match()
                self.display.tile('column')


            images = [corrected]
            labels = ['science']
            for i in np.arange(len(labels)):
                self.set(exposureprefix, labels[i])
                self.data = images[i]

            self.speak('total corrected flux is {0}'.format(np.sum(image - corrected)))

            lostflux = np.sum(image - corrected, axis=0)
            try:
                self.cosmicplot.set_ydata(lostflux)
            except AttributeError:
                plt.figure('cosmic ray rejection', figsize=(5, 3), dpi=100)
                self.axcr = plt.subplot()
                self.axcr.cla()
                self.cosmicplot = self.axcr.plot(lostflux, color='Sienna')[0]
                self.axcr.set_ylim(1.0, 1e8)
                self.axcr.set_xlim(-1, len(lostflux)+1)
                self.axcr.set_yscale('log')
                self.axcr.set_ylabel('Total Flux Rejected (e-/column)')
                self.axcr.set_xlabel('Column (pixels)')
                plt.tight_layout()

            self.axcr.set_title('Science {0}'.format(exposureprefix))
            plt.draw()
            np.save(self.cosmicsFilename, lostflux)
            plt.savefig(self.cosmicsFilename.replace('.npy', '.png'))
            self.speak('saved cosmic ray rejection checks to {0}'.format(cosmics_directory))
            #self.input('thoughts on CR?')

    def loadImages(self, exposureprefixes, imageType=None):
        '''Load a series of CCD images, returning them as a cube.'''

        # if n is a single element array, just return one image
        try:
            exposureprefixes[1]
        except (ValueError, IndexError):
            self.set(exposureprefixes[0])
            return self.readData(imageType=imageType)

        # loop over the image numbers, and read them
        images = []
        for i in range(len(exposureprefixes)):
            images.append(self.readData(exposureprefixes[i], imageType))

        # convert a list of images into an array
        cube = np.array(images)

        # return that array
        return cube

    def createStitched(self):
        '''Create and load a stitched CCD image, given a (stored) file prefix.'''

        # print status
        self.speak("creating a stitched image for {0}".format(self.stitched_filename))

        # provide different options for different kinds of images
        if self.imageType == 'bias':
            self.flags['subtractbias'] = False
            self.flags['subtractdark'] = False
            self.flags['multiplygain'] = False
        elif self.imageType == 'dark':
            self.flags['subtractbias'] = True
            self.flags['subtractdark'] = False
            self.flags['multiplygain'] = False
        elif self.imageType == 'FlatInADU':
            self.flags['subtractbias'] = True
            self.flags['subtractdark'] = True
            self.flags['multiplygain'] = False
        else:
            self.flags['subtractbias'] = True
            self.flags['subtractdark'] = True
            self.flags['multiplygain'] = True

        # don't restitch if unnecessary
        if os.path.exists(self.stitched_filename):
            self.speak("{0} has already been stitched".format(self.name))
        else:
            # process the two halves separately, and then smush them together
            filenames = [os.path.join(self.obs.night.dataDirectory, f) for f in self.instrument.prefix2files(self.exposureprefix)]

            # load the (only) image
            stitched, header = self.instrument.loadSingleCCD(filenames)

            if self.visualize:
                tempstitched = stitched

            if self.visualize:
                self.display.one(stitched, clobber=True)
                self.input('This is the raw stitched image; press enter to continue.')

            # subtract bias
            if self.flags['subtractbias']:
                self.speak("subtracting bias image")

                # pull the bias from the calibraion object
                stitched -= self.calib.bias()

            if self.visualize:
                self.display.one(stitched, clobber=True)
                self.input('after subtracting bias')

            # normalize darks by exposure time
            if self.imageType == 'dark':
                # if we're creating a dark, normalize to e/s
                stitched /= self.instrument.darkexptime(header)

            # subtract dark
            if self.flags['subtractdark']:
                self.speak("subtracting dark image")
                # pull the dark image from the calibration object; multiply by the effective exposure time for the image we're trying to calibrate
                stitched -= self.calib.dark()*self.instrument.darkexptime(header)

            if self.visualize:
                self.display.one(stitched, clobber=True)
                self.visualize = self.input('after subtracting dark; type [s] to stop showing these').lower() != 's'

            # divide by the gain (pulled from the header)
            if self.flags['multiplygain']:

                # if there are multiple gains for different parts of the chip, you may need to make gain() return an image of the appropriate shape
                gain = self.instrument.gain(header)
                self.speak("multiplying by gain of {0} e-/ADU".format(np.unique(np.asarray(gain))))
                # convert from ADU to electrons
                stitched *= gain

            if self.visualize:
                self.display.one(stitched, clobber=True)
                self.visualize = self.input('after multiplying by gain; type [s] to stop showing these').lower() != 's'

            # put the stitched image into the CCD's memory
            self.data = stitched

            # find and reject cosmics based on nearby images in time
            if self.instrument.zapcosmics:
                if self.imageType == 'science':
                    self.rejectCosmicRays() # KLUDGE -- I'm pretty sure this shouldn't be used

            # write out the image to a stitched image
            writeFitsData(self.data, self.stitched_filename)
            self.speak("stitched and saved {0}".format(self.name))
            assert(np.isfinite(self.data).any())
