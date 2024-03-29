from .Spectrograph import *

class LDSS3C(Spectrograph):
    '''
    This is an LDSS3C instrument. It contains all instrument specific
    parameters and processes. Some procedures are inherited generally
    from Spectrograph, so functions that you expect to be common to
    many instruments should likely be coded there.
    '''

    # a string for the instrument name
    name = 'LDSS3C'

    # are the slits in a "mask" (with different locations for every star)
    #               or a "longslit" (with one location for each star?)
    slitstyle = 'mask'

    # file search patten to get a *single* fileprefix for each exposure
    # LDSS3C has two amplifier readouts ('c1.fits' and 'c2.fits');
    #  this pulls the fileprefix just from c1
    basicpattern = 'ccd*c1.fits'

    # which header of the fits file contains the header with useful information
    fitsextensionforheader = 0

    # what keys should be included in the nightly logs for this instrument?
    # this is everything you might want to have access to at some point
    keysforlogheader = [    'ut-date',
                            'ut-time',
                            'object',
                            'exptime',
                            'aperture',
                            'filter',
                            'grism',
                            'airmass',
                            'filename',
                            'ra-d', 'dec-d',
                            'exptype',
                            'binning',
                            'speed',
                            'gain',
                            'comment']

    # what keys should make it into condensed summary logs?
    # these are things you'll want for deciding what type each file is
    keysforsummary = [      'fileprefix',
                            'object',
                            'exptime',
                            'aperture',
                            'filter',
                            'grism',
                            'airmass']

    # what keys do we want to store associated with a science timeseries?
    # these will show up, ultimately, in the 'temporal' key of a cube
    keysfortimeseries = [   'date-obs',
                            'ut-date',
                            'ut-time',
                            'ut-end',
                            'scale',
                            'gain',
                            'epoch',
                            'airmass',
                            'ha',
                            'exptime',
                            'tempccd',
                            'templdss',
                            'focus',
                            'rotangle',
                            'rotatore']

    globallinekeys = ['airmass', 'rotatore']
    # these keys are useful to search for guessing the filetype
    # for LDSS3C, the usful information is in the "object" key of the header
    # for other instruments, I could imagine "comments" being useful sometimes
    keytosearch = 'object'

    # by what key should files be sorted in the summaries?
    summarysortkey = 'fileprefix'

    # within that header key, what words do we search for?
    wordstosearchfor = { 'dark':['dark', ''],
                         'bias':['bias', 'Bias', ''],
                         'flat':['quartz', 'flat', 'lamp', ''],
                           'He':['He', 'helium'],
                           'Ne':['Ne', 'neon'],
                           'Ar':['Ar', 'argon'],
                      'science':['sciece', 'spec']}

    wordstoavoid  =    { 'dark':[],
                         'bias':[],
                         'flat':[],
                           'He':[],
                           'Ne':[],
                           'Ar':['dark', 'quartz'],
                      'science':['science']}

    def __repr__(self):
        '''
        How should this object be represented as a string?
        '''
        return '<Spectrograph {}>'.format(self.name)

    def __init__(self, grism='vph-red'):

        # what's the name of this instrument?
        self.name = 'LDSS3C'

        # where is it located? (needed for BJD calculation)
        self.telescope = 'Magellan'
        self.sitename = 'LCO'

        #try:
        #self.observatory = coord.EarthLocation.of_site(self.sitename)
        # what used to be here:
        #self.observatory = coord.EarthLocation.from_geodetic(-4.71333*u.hourangle, -29.00833*u.deg, 2282.0*u.m)
        # testing with revised location from LCO.CL website
        self.observatory = coord.EarthLocation.from_geodetic(lon=coord.Angle('-70° 41′ 33.36″'), lat=coord.Latitude('-29° 0′ 52.56″'), height=2380*u.m)
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
        #self.gains = np.array([1.67, 1.45])  # gain from 2017, fast read out (user manual)
        #self.gains = np.array([1.570, 1.383])  # new gain from working on instrument, fast read out (20 June 2019)
        #self.gains = np.array([2.95, 2.56])  # gain from 2019, turbo read out mode (August 2019)
        self.gains = np.array([2.84, 2.52])  # gain from 2021, turbo read out mode (March 2021)
        self.binning = 1

        # what area of the detector contains real data? (for individual amplifiers
        self.dataleft = 0
        if   self.binning == 1: self.dataright = 1024
        elif self.binning == 2: self.dataright = 512
        self.databottom = 0
        if   self.binning == 1: self.datatop = 4096
        elif self.binning == 2: self.datatop = 2048

        # set up the size of the image
        self.xsize = self.namps*(self.dataright - self.dataleft)
        self.ysize = (self.datatop - self.databottom)

        # what are the calibrations we should expect
        self.detectorcalibrations = ['dark', 'bias', 'flat']

        # how many stitched images can we hold in memory?
        self.maximumimagesinmemory = 128
        self.maximumimagesinmemoryforscience = 75


    def setupDisperser(self):
        '''
        Setup the basics for the disperser.
        '''

        # what is the name and type of the disperser
        self.dispersertype = 'grism'
        self.disperser = self.grism

        # define a uniform grid of wavelengths for supersampling onto, later
        if self.grism == 'vph-all':
            self.uniformwavelengths = np.arange(5000, 9800)
            self.alignmentranges = dict(    #Hbeta=(4750,5050),
                                            Halpha=(6425,6725),
                                            O2B=(6850,6950),
                                            O2A=(7550,7750),
                                            #Ca_triplet=(8450,8725),
                                            H2O=(9300, 9500)
                                            )
        elif self.grism == 'vph-red':
            self.uniformwavelengths = np.arange(6000, 10500)
            self.alignmentranges = dict(    #Hbeta=(4750,5050),
                                            Halpha=(6425,6725),
                                            O2B=(6850,6950),
                                            O2A=(7550,7750),
                                            Ca_triplet=(8450,8725),
                                            H2O=(9300, 9500)
                                            )

        # the available arc lamps for wavelength calibration
        self.arclamps = ['He', 'Ne', 'Ar']

        # set up the wavelength calibration paths and files
        self.disperserDataDirectory = os.path.join(mosasaurusdirectory,
                                                'data/',
                                                self.name + '/',
                                                self.disperser + '/')
        self.wavelength2pixelsFile = os.path.join(self.disperserDataDirectory,
                '{0}_wavelength_identifications.txt'.format(self.grism))

        self.wavelengthsFile = os.path.join(self.disperserDataDirectory,
                'HeNeAr.txt')

        if self.binning == 2:
            self.offsetBetweenReferenceAndWavelengthIDs = -1024
        elif self.binning == 1: 
            self.offsetBetweenReferenceAndWavelengthIDs = -2048
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

    def gain(self, header):

        zeros = [np.zeros((self.datatop, self.dataright)) for i in self.gains]
        gains = [zeros[i] + self.gains[i] for i in range(len(zeros))]
        gainimage = self.stitchChips(gains)
        return gainimage

    def exptime(self, header):
        return header[0]['EXPTIME']

    def darkexptime(self, header):
        return header[0]['EXPTIME']

    def stitchChips(self, listOfChips):
        # for now just working with chip8
        return np.hstack([listOfChips[0], np.fliplr(listOfChips[1])])

    def loadOverscanTrimHalfCCD(self, filename):
        '''Open one half of an amplifier of a chip, subtract the overscan, and trim.'''

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

    def loadSingleCCD(self, filenames):
        '''
        Load an IMACS image; subtract and trim its overscan.
        In general, this function should load and return a single image.
        If the detector uses multiple amplifiers to read out different parts of the same chip,
        this function should stitch those sections together.

        Parameters
        ----------
        filenames: list
            a list of relevant filenames (e.g. multiple amplifiers)

        Returns
        -------
        image: array
            a overscan-trimmed CCD image
        '''

        # load the chips
        c_ = np.array([self.loadOverscanTrimHalfCCD(f) for f in filenames])
        c_data = c_[:,0]    # list of chip data
        c_header = c_[:,1]  # list of chip header

        # stitch the CCD's together
        stitched = self.stitchChips(c_data)

        return stitched, c_header
