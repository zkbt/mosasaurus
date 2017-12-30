from .Spectrograph import *

class LDSS3C(Spectrograph):

    name = 'LDSS3C'

    basicpattern = 'ccd*c1.fits'

    # which header of the fits file contains the header with useful information
    fitsextensionforheader = 0

    # what keys should be included in the nightly logs for this instrument?
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
    keysforsummary = [      'fileprefix',
                            'object',
                            'exptime',
                            'aperture',
                            'filter',
                            'grism',
                            'airmass']

    # what keys do we want to store associated with a science timeseries?
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


    # these keys are useful to search for guessing the filetype
    keytosearch = 'object'

    # within those keys, what words do we search for?
    wordstosearchfor = {'dark':['dark'],
                         'bias':['bias'],
                         'flat':['quartz', 'flat'],
                         'He':['He', 'helium'],
                         'Ne':['Ne', 'neon'],
                         'Ar':['Ar', 'argon']}

    def __repr__(self):
        '''How should this object be represented as a string?'''
        return '<Spectrograph {}>'.format(self.name)

    def findDarks(self, night):
        '''Identify the dark exposures.'''
        match = night.find( wordstolookfor = self.wordstosearchfor['dark'],
                            placetolook = self.keytosearch)
        return match

    def findBiases(self, night):
        '''Identify the bias exposures.'''
        match = night.find( wordstolookfor = self.wordstosearchfor['bias'],
                            placetolook = self.keytosearch)
        return match

    def findHe(self, night):
        '''Identify the He exposures.'''
        match = night.find( wordstolookfor = self.wordstosearchfor['He'],
                            placetolook = self.keytosearch)
        return match

    def findNe(self, night):
        '''Identify the Ne exposures.'''
        match = night.find( wordstolookfor = self.wordstosearchfor['Ne'],
                            placetolook = self.keytosearch)
        return match

    def findAr(self, night):
        '''Identify the Ar exposures.'''
        match = night.find( wordstolookfor = self.wordstosearchfor['Ar'],
                            placetolook = self.keytosearch)
        return match

    def __init__(self, grism='vph-red'):

        # what's the name of this instrument?
        self.name = 'LDSS3C'

        # where is it located? (needed for BJD calculation)
        self.telescope = 'Magellan'
        self.sitename = 'LCO'
        self.observatory = coord.EarthLocation.of_site(self.sitename)

        '''
        # removed this once I realized astropy knew about LCO
        self.observatory = dict(
            name="Las Campanas Observatory",
            timezone="Chilean",
            standardzone = 4.0*astropy.units.hour,
            usedaylightsaving = -1,
            longitude = 4.71333*astropy.units.hourangle,
            latitude = -29.00833*astropy.units.deg,
            elevsea = 2282.0*astropy.units.m,
            elev = 2282.0*astropy.units.m, # /* for ocean horizon, not Andes! */
            )
        '''

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

        # what are the calibrations we should expect
        self.detectorcalibrations = ['dark', 'bias', 'flat']

    def setupDisperser(self):
        '''
        Setup the basics for the disperser.
        '''

        # what is the name and type of the disperser
        self.dispersertype = 'grism'
        self.disperser = self.grism

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

    def setupExtraction(self):
        '''
        Setup the default extraction parameters associated with this instrument.
        '''
        self.extractiondefaults = {}

        # the geometry of the pixels to consider for each extraction
        # how many pixels in the spatial direction should analysis extend?
        self.extractiondefaults['spatialsubarray'] = 50
        # how far (in pixels) does spectrum extend away from direct image position
        self.extractiondefaults['wavelengthredward'] = np.inf
        self.extractiondefaults['wavelengthblueward'] = np.inf


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
        zachopy.utils.mkdir(self.dataDirectory)

        # working directory is where in-progress files + results are created
        self.workingDirectory = os.path.join(self.baseDirectory,
                                                workingDirectory)
        zachopy.utils.mkdir(self.workingDirectory)

        # extraction directory is where extracted spectra land
        self.extractionDirectory = os.path.join(self.baseDirectory,
                                                extractionDirectory)
        zachopy.utils.mkdir(self.extractionDirectory)


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

    def fileprefix(self, n):
        '''
        Feed in a CCD number,
        spit out the file prefix for
        that CCD amplifier pair.
        '''
        try:
          return [self.dataDirectory + 'ccd{0:04}'.format(x) for x in n]
        except TypeError:
          return self.dataDirectory + 'ccd{0:04}'.format(n)

    def file2prefix(self, filename):
        '''
        This function returns a shortened fileprefix from a given filename.
        '''
        tail = os.path.split(filename)[-1]

        # LDSS3C is split into two amplifiers, let's pull out the prefix
        return tail.replace('c2.fits', '').replace('c1.fits', '')

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

#def identifyImageNumbers(self, lookingfor)
