from .Spectrograph import *

class IMACS(Spectrograph):

    name = 'IMACS'

    basicpattern = 'ift*c8.fits' # for now we will just extract chip 8; stitching is hard

    # which header of the fits file contains the header with useful information
    fitsextensionforheader = 0

    # what keys should be included in the nightly logs for this instrument?
    keysforlogheader = [    'ut-date',
                            'ut-time',
                            'object',
                            'exptime',
                            'slitmask',
                            'filter',
                            'dispersr',
                            'airmass',
                            'filename',
                            'ra-d', 'dec-d',
                            'exptype',
                            'binning',
                            'speed',
                            'ccdgain',
                            'comment']

    # what keys should make it into condensed summary logs?
    keysforsummary = [      'fileprefix',
                            'object',
                            'exptime',
                            'slitmask',
                            'filter',
                            'dispersr',
                            'airmass']

    # what keys do we want to store associated with a science timeseries?
    keysfortimeseries = [   'date-obs',
                            'ut-date',
                            'ut-time',
                            'ut-end',
                            'scale',
                            'ccdgain',
                            'epoch',
                            'airmass',
                            'ha-str',
                            'exptime',
                            'tempccd8',
                            'tempstr',
                            'detfocus',
                            'rotangle',
                            'rotatore']


    # these keys are useful to search for guessing the filetype
    keytosearch = 'object'

    # by what key should files be sorted in the summaries?
    summarysortkey = 'fileprefix'

    # within those keys, what words do we search for?
    wordstosearchfor = {'dark':['dark'],
                        'bias':['bias'],
                        'flat':['quartz', 'flat'],
                        'He':['He', 'helium'],
                        'Ne':['Ne', 'neon'],
                        'Ar':['Ar', 'argon']}

    wordstoavoid  =    { 'dark':[],
                         'bias':[],
                         'flat':[],
                           'He':[],
                           'Ne':[],
                           'Ar':['dark', 'quartz']}

    def __repr__(self):
        '''How should this object be represented as a string?'''
        return '<Spectrograph {}>'.format(self.name)

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
        self.name = 'IMACS'

        # where is it located? (needed for BJD calculation)
        self.telescope = 'Magellan'
        self.sitename = 'LCO'

        #try:
        #self.observatory = coord.EarthLocation.of_site(self.sitename)
        self.observatory = coord.EarthLocation.from_geodetic(lon=coord.Angle('-70° 41′ 33.36″'), lat=coord.Latitude('-29° 0′ 52.56″'), height=2380*u.m)
        #EarthLocation(1845655.49905341*m, -5270856.2947176*m, -3075330.77760682*m)

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

        # what grism is being used ['gri-300-26.7']
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
        #self.namps = 1  # for now we are only extracting chip 8
        #self.gains = np.array([1.50])
        self.namps = 8
        self.gains = np.array([1.47, 1.47, 1.58, 1.50, 1.52, 1.54, 1.48, 1.50])
        #self.namps = 4
        #self.gains = np.array([1.52, 1.54, 1.48, 1.50])

        self.binning = 2

        # what area of the detector contains real data? (for individual amplifiers
        self.dataleft = 0
        self.dataright = 1024
        self.databottom = 0
        self.datatop = 2048

        # set up the size of the image
        self.xsize = 4096
        self.ysize = 4096

        # what are the calibrations we should expect
        self.detectorcalibrations = ['flat'] # for some data sets where I didn't take biases or flats (you shouldn't need them for IMACS)
        #self.detectorcalibrations = ['dark', 'bias', 'flat']

        # how many stitched images can we hold in memory?
        self.maximumimagesinmemory = 75
        self.maximumimagesinmemoryforscience = 32

    def setupDisperser(self):
        '''
        Setup the basics for the disperser.
        '''

        # what is the name and type of the disperser
        self.dispersertype = 'grism'
        self.disperser = self.grism

        # define a uniform grid of wavelengths for supersampling onto, later
        if self.grism == 'gri-150-10.8':
            self.uniformwavelengths = np.arange(5000, 10000)
            self.alignmentranges = {    r'$O_2$ - A':(7500,7800),
                                        r'Ca triplet':(8450,8750),
                                        r'$H_2O$':(9200, 9500),
                                            }

        if self.grism == 'gri-300-26.7':
            self.uniformwavelengths = np.arange(5000, 10000)
            self.alignmentranges = {    r'$O_2$ - A':(7500,7800),
                                        r'Ca triplet':(8450,8750),
                                        r'$H_2O$':(9200, 9500),
                                            }
            
        self.offsetBetweenReferenceAndWavelengthIDs = 0.

        # the available arc lamps for wavelength calibration
        self.arclamps = ['He', 'Ne', 'Ar']

        # set up the wavelength calibration paths and files
        self.disperserDirectory = os.path.join(mosasaurusdirectory,
                                                'data/',
                                                self.name + '/',
                                                self.disperser + '/')
        self.wavelength2pixelsFile = os.path.join(self.disperserDirectory,
                '{0}_wavelength_identifications.txt'.format(self.grism))    # at this stage this .txt file is NOT the correct one for IMACS

        self.wavelengthsFile = os.path.join(self.disperserDirectory,
                'HeNeAr.txt')

        # offset to where wavelengths were idenified

        # find the peak of the combined correlation function
        #if self.aperture.obs.instrument == 'LDSS3C':
        #    self.peakoffset = -1024 # KLUDGE KLUDGE KLUDGE! np.where(self.corre['combined'] == self.corre['combined'].max())[0][0] - len(x)
        #    # (old?) to convert: len(x) - xPeak = x + peakoffset
        #elif self.aperture.obs.instrument == 'IMACS': self.peakoffset = -75  # also a kludge

    def setupExtraction(self):
        '''
        Setup the default extraction parameters associated with this instrument.
        '''
        self.extractiondefaults = {}

        # the geometry of the pixels to consider for each extraction
        # how many pixels in the spatial direction should analysis extend?
        self.extractiondefaults['spatialsubarray'] = 50
        # how far (in pixels) does spectrum extend away from direct image position
        if self.grism == 'gri-150-10.8':
            self.extractiondefaults['stampwavelengthredward'] = 1000
            self.extractiondefaults['stampwavelengthblueward'] = 500

        if self.grism == 'gri-300-26.7':
            self.extractiondefaults['stampwavelengthredward'] = np.inf
            self.extractiondefaults['stampwavelengthblueward'] = 1200



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

        # we're going to cross-correlate the spectra within this wavelength range
        self.extractiondefaults['correlationAnchors'] = [8498.0, 8542.0, 8662.0]
        self.extractiondefaults['correlationRange'] = [8350, 8800]
        self.extractiondefaults['correlationSmooth'] = 2


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

        return tail.replace('c1.fits', '').replace('c2.fits', '').replace('c3.fits', '').replace('c4.fits', '').replace('c5.fits', '').replace('c6.fits', '').replace('c7.fits', '').replace('c8.fits', '')

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
        return [prefix + 'c1.fits', prefix + 'c2.fits', prefix + 'c3.fits', prefix + 'c4.fits', prefix + 'c5.fits', prefix + 'c6.fits', prefix + 'c7.fits', prefix + 'c8.fits']

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

        bottomrow = np.hstack((np.flipud(listOfChips[6]), np.flipud(listOfChips[7]), np.flipud(listOfChips[4]), np.flipud(listOfChips[5])))
        toprow = np.hstack((np.fliplr(listOfChips[3]), np.fliplr(listOfChips[2]), np.fliplr(listOfChips[1]), np.fliplr(listOfChips[0])))
        return np.vstack((toprow, bottomrow))

    def loadOverscanTrimHalfCCD(self, filename):
        '''Open one half of an amplifier of a chip, subtract the overscan, and trim.'''
        ### FIX ME ### -- this should all be moved into the Spectrograph definition

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


