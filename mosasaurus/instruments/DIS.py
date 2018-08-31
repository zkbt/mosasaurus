from .Spectrograph import *

class DIS(Spectrograph):
    '''
    This is an DIS instrument. It contains all instrument specific
    parameters and processes. Some procedures are inherited generally
    from Spectrograph, so functions that you expect to be common to
    many instruments should likely be coded there.
    '''

    # a string for the instrument name
    name = 'DIS'

    # file search patten to get a *single* fileprefix for each exposure
    # DIS has two camera readouts ('*b.fits' and '*r.fits' for blue and red);
    # this pulls the fileprefix just from blue
    basicpattern = '*b.fits'

    # which header of the fits file contains the header with useful information
    fitsextensionforheader = 0

    # what keys should be included in the nightly logs for this instrument?
    # this is everything you might want to have access to at some point
    keysforlogheader = [    'DATE-OBS',
                            'OBJNAME',
                            'EXPTIME',
                            'SLITMASK',
                            'FILTER',
                            'GRATING',
                            'AIRMASS',
                            'FILENAME',
                            'RA', 'DEC',
                            'IMAGETYP',
                            'DETECTOR',
                            'CCDBIN1', 'CCDBIN2',
                            'GAIN', 'RDNOISE']

    # what keys should make it into condensed summary logs?
    # these are things you'll want for deciding what type each file is
    keysforsummary = [      'FILENAME',
                            'OBJNAME',
                            'EXPTIME',
                            'SLITMASK',
                            'GRATING',
                            'AIRMASS']

    # what keys do we want to store associated with a science timeseries?
    # these will show up, ultimately, in the 'temporal' key of a cube
    keysfortimeseries = [   'DATE-OBS',
                            'LST',
                            'OBJANGLE',
                            'AIRPRESS',
                            'HUMIDITY',
                            'TELAZ', 'TELALT', 'TELROT', 'TELFOCUS',
                            'AIRMASS',
                            'EXPTIME',
                            'DARKTIME',
                            'CCDTEMP', 'CCDHEAT']


    # these keys are useful to search for guessing the filetype
    # for LDSS3C, the usful information is in the "object" key of the header
    # for other instruments, I could imagine "comments" being useful sometimes
    keytosearch = 'IMAGETYP'

    # within that header key, what words do we search for?
    wordstosearchfor = { 'dark':['dark'],
                         'bias':['bias'],
                         'flat':['quartz', 'flat'],
                           'He':['He', 'helium'],
                           'Ne':['Ne', 'neon'],
                           'Ar':['Ar', 'argon']}



    def __init__(self, grism='vph-red'):

        # what's the name of this instrument?
        self.name = 'LDSS3C'

        # where is it located? (needed for BJD calculation)
        self.telescope = 'Magellan'
        self.sitename = 'LCO'

        #try:
        #self.observatory = coord.EarthLocation.of_site(self.sitename)
        self.observatory = coord.EarthLocation.from_geodetic(-4.71333*u.hourangle, -29.00833*u.deg, 2282.0*u.m)
        #EarthLocation(1845655.49905341*m, -5270856.2947176*m, -3075330.77760682*m)


        # what grism is being used ['vph-red', 'vph-all', 'vph-blue']
        self.grism = grism.lower()

        # run the setup scripts, once these basics are defined
        Spectrograph.__init__(self)

    def setupDetector(self):
        '''
        Setup the basics for the detector.
        (These can be defaults can by modified by later changing
        the attributes, like for example:

            l = LDSS3C()
            l.gains = np.array([42.0, 54.0])
        '''

        # basic information about the amplifiers
        self.namps = 2
        self.gains = np.array([1.72, 1.49])
        self.binning = 2

        # what area of the detector contains real data? (for individual amplifiers
        self.dataleft = 0
        self.dataright = 512
        self.databottom = 0
        self.datatop = 2048

        # set up the size of the image
        self.xsize = self.namps*(self.dataright - self.dataleft)
        self.ysize = (self.datatop - self.databottom)

        # what are the calibrations we should expect
        self.detectorcalibrations = ['dark', 'bias', 'flat']

    def setupDisperser(self):
        '''
        Setup the basics for the disperser.
        '''

        # what is the name and type of the disperser
        self.dispersertype = 'grism'
        self.disperser = self.grism

        # define a uniform grid of wavelengths for supersampling onto, later
        if self.grism == 'vph-all':
            self.uniformwavelengths = np.arange(4000, 10500)
            self.alignmentranges = {    r'$H\beta$':(4750,5050),
                                        r'$H\alpha$':(6425,6725),
                                        r'$O_2$ - B':(6750,7050),
                                        r'$O_2$ - A':(7500,7800),
                                        r'Ca triplet':(8450,8750),
                                        r'$H_2O$':(9200, 9700),
                                            }
        elif self.grism == 'vph-red':
            self.uniformwavelengths = np.arange(6000, 10500)
            self.alignmentranges = dict(    O2=(7580, 7650),
                                            Ca1=(8490, 8525),
                                            Ca2=(8535, 8580),
                                            Ca3=(8650, 8700),
                                            H2O=(9300, 9700)
                                            )


        # the available arc lamps for wavelength calibration
        self.arclamps = ['He', 'Ne', 'Ar']

        # set up the wavelength calibration paths and files
        self.disperserDirectory = os.path.join(mosasaurusdirectory,
                                                'data/',
                                                self.name + '/',
                                                self.disperser + '/')
        self.wavelength2pixelsFile = os.path.join(self.disperserDirectory,
                '{0}_wavelength_identifications.txt'.format(self.grism))

        self.wavelengthsFile = os.path.join(self.disperserDirectory,
                'HeNeAr.txt')

        if self.binning == 2:
            self.offsetBetweenReferenceAndWavelengthIDs = -1024

        # find the peak of the combined correlation function
        #if self.aperture.obs.instrument == 'LDSS3C':
        #    self.offsetBetweenReferenceAndWavelengthIDs = -1024 # KLUDGE KLUDGE KLUDGE! np.where(self.corre['combined'] == self.corre['combined'].max())[0][0] - len(x)
        #    # (old?) to convert: len(x) - xPeak = x + offsetBetweenReferenceAndWavelengthIDs
        #elif self.aperture.obs.instrument == 'IMACS': self.offsetBetweenReferenceAndWavelengthIDs = -75  # also a kludge

    def setupExtraction(self):
        '''
        Setup the default extraction parameters associated with this instrument.
        '''
        self.extractiondefaults = {}

        # the geometry of the pixels to consider for each extraction
        # how many pixels in the spatial direction should analysis extend?
        self.extractiondefaults['spatialsubarray'] = 50
        # how far (in pixels) does spectrum extend away from direct image position
        self.extractiondefaults['stampwavelengthredward'] = np.inf
        self.extractiondefaults['stampwavelengthblueward'] = np.inf


        # setup the default initial extraction geometry
        #  (some of these these can be modified interactively later)
        # what are the aperture widths to consider?
        self.extractiondefaults['narrowest'] = 2
        self.extractiondefaults['widest'] = 12
        self.extractiondefaults['numberofapertures'] = 6
        # what order polynomial to use for trace shape?
        self.extractiondefaults['traceOrder'] = 2
        # initial guess for the width of the sky regions
        self.extractiondefaults['skyWidth'] = 10
        # required minimum gap between extraction and sky apertures
        self.extractiondefaults['skyGap'] = 2
        # should we try to zap cosmic rays?
        self.extractiondefaults['zapcosmics'] = False

        # what are the kinds of images extractions can work with
        self.extractables = ['science', 'reference']


    def setupDirectories(self,
            baseDirectory='/Users/zkbt/Cosmos/Data/Magellan/LDSS3/',
            dataDirectory='data/',
            workingDirectory='working/',
            extractionDirectory='extraction/'):
        '''
        Setup the basic file structure associated with this instrument.
        If you store your data anywhere other than the defaults listed here,
        then adjust the directory structure with:
            i = LDSS3C()
            i.setupDirectories("/absolute/path/to/base/directory/")
        '''

        # absolute path to where the data are located
        self.baseDirectory = baseDirectory

        # data directory is where raw (not-to-be-modified) files are stored
        self.dataDirectory = os.path.join(self.baseDirectory,
                                                dataDirectory)
        mkdir(self.dataDirectory)

        # working directory is where in-progress files + results are created
        self.workingDirectory = os.path.join(self.baseDirectory,
                                                workingDirectory)
        mkdir(self.workingDirectory)

        # extraction directory is where extracted spectra land
        self.extractionDirectory = os.path.join(self.baseDirectory,
                                                extractionDirectory)
        mkdir(self.extractionDirectory)


    def extractMidExposureTimes(self, headers):
        '''
        For an astropy table of extracted header values,
        extract the mid-exposure times in JD_UTC.
        '''

        # stitch together date+time strings
        timestrings = ['{0} {1}'.format(row['ut-date'], row['ut-time']) for row in headers]

        # calculate a JD from these times (and adding half the exposure time, assuming ut-time is the start of the exposure)
        starttimes = astropy.time.Time(timestrings, format='iso', scale='utc', location=self.observatory)

        # mid-exposure
        times_earth = starttimes + 0.5*headers['exptime']*u.second

        # return as astropy times, in the UTC system, at the location of the telescope
        return times_earth


    def file2prefix(self, filename):
        '''
        This function returns a shortened fileprefix from a given filename.
        '''
        tail = os.path.split(filename)[-1]

        # LDSS3C is split into two amplifiers, let's pull out the prefix
        return tail.replace('c2.fits', '').replace('c1.fits', '')

    def prefix2number(self, prefix):
        '''
        This function returns a CCD number (not necessarily starting from 0)
        from a fileprefix.
        '''
        return np.int(prefix[-4:])

    def prefix2files(self, prefix):
        '''
        This function returns a list of filenames (without complete path)
        that are associated with this given prefix.
        '''
        return [prefix + 'c1.fits', prefix + 'c2.fits']


    def loadOverscanTrimHalfCCD(self, filename):
        '''
        Open one half of an LDSS3 CCD, subtract the overscan, and trim.
        '''

        # open the FITS file, split into header and data
        hdu = astropy.io.fits.open(filename)
        header = hdu[0].header
        data = readFitsData(filename)

        # take the parts of CCD exposed to light
        goodData = data[self.databottom:self.datatop,self.dataleft:self.dataright]
        goodBias = data[self.databottom:self.datatop,self.dataright:]

        # estimate the 1D bias (and drawdown, etc...) from the overscan
        biasEstimate = np.median(goodBias, axis=1)
        biasImage = np.ones(goodData.shape)*biasEstimate[:,np.newaxis]

        return (goodData - biasImage), header

    def createStitched(self, ccd):
        '''Create and load a stitched CCD image, given a file prefix.'''

        # print status
        self.speak("creating a stitched image for {0}".format(ccd.stitched_filename))

        # provide different options for different kinds of images
        if ccd.imageType == 'bias':
            ccd.flags['subtractbias'] = False
            ccd.flags['subtractdark'] = False
            ccd.flags['multiplygain'] = False
            ccd.flags['subtractcrosstalk'] = False
        elif ccd.imageType == 'dark':
            ccd.flags['subtractbias'] = True
            ccd.flags['subtractdark'] = False
            ccd.flags['multiplygain'] = False
            ccd.flags['subtractcrosstalk'] = False
        elif ccd.imageType == 'FlatInADU':
            ccd.flags['subtractbias'] = True
            ccd.flags['subtractdark'] = True
            ccd.flags['multiplygain'] = False
            ccd.flags['subtractcrosstalk'] = False
        else:
            ccd.flags['subtractbias'] = True
            ccd.flags['subtractdark'] = True
            ccd.flags['multiplygain'] = True
            ccd.flags['subtractcrosstalk'] = True

        # don't restitch if unnecessary
        if os.path.exists(ccd.stitched_filename):
            self.speak("{0} has already been stitched".format(self.name))
        else:
            # process the two halves separately, and then smush them together
            filenames = [os.path.join(self.obs.night.dataDirectory, f) for f in self.prefix2files(ccd.exposureprefix)]


            # load the two halves
            c1data, c1header = self.loadOverscanTrimHalfCCD(filenames[0])
            c2data, c2header = self.loadOverscanTrimHalfCCD(filenames[1])

            if ccd.visualize:
                tempstitched = np.hstack([c1data, np.fliplr(c2data)])

            if ccd.flags['subtractcrosstalk']:
                    # is this possible?
                    pass

            # stitch the CCD's together
            stitched = np.hstack([c1data, np.fliplr(c2data)])

            if ccd.visualize:
                ccd.display.one(stitched, clobber=True)
                self.input('This is the raw stitched image; press enter to continue.')

            # subtract bias
            if ccd.flags['subtractbias']:
                self.speak("subtracting bias image")
                stitched -= ccd.calib.bias()

            if ccd.visualize:
                ccd.display.one(stitched, clobber=True)
                self.input('after subtracting bias')

            # normalize darks by exposure time
            if ccd.imageType == 'dark':
                stitched /= c1header['EXPTIME']

            # subtract dark
            if ccd.flags['subtractdark']:
                self.speak("subtracting dark image")
                stitched -= ccd.calib.dark()*c1header['EXPTIME']

            if ccd.visualize:
                ccd.display.one(stitched, clobber=True)
                ccd.visualize = self.input('after subtracting dark; type [s] to stop showing these').lower() != 's'

            # divide by the gain (KLUDGE! make sure these are the best estimates!)
            if ccd.flags['multiplygain']:

                try:
                    self.gains
                except AttributeError:
                        self.ccd.estimateGain()

                self.speak("multiplying by gains of {0} e-/ADU".format(self.gains))
                gain1 = np.zeros_like(c1data) + self.gains[0]
                gain2 = np.zeros_like(c2data)+ self.gains[1]
                gainimage = np.hstack([gain1, np.fliplr(gain2)])
                stitched *= gainimage

            if ccd.visualize:
                ccd.display.one(stitched, clobber=True)
                ccd.visualize = self.input('after multiplying by gain; type [s] to stop showing these').lower() != 's'

            # put the stitched image into the CCD's memory
            ccd.data = stitched

            # find and reject cosmics based on nearby images in time
            if self.zapcosmics:
                if ccd.imageType == 'science':
                    ccd.rejectCosmicRays() # KLUDGE -- I'm pretty sure this shouldn't be used

            # write out the image to a stitched image
            writeFitsData(ccd.data, ccd.stitched_filename)
            self.speak("stitched and saved {0}".format(ccd.name))


#def identifyImageNumbers(self, lookingfor)
