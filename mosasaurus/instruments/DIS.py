from .Spectrograph import *

class DIS(Spectrograph):
    '''
    This is an DIS instrument. It contains all instrument specific
    parameters and processes. Some procedures are inherited generally
    from Spectrograph, so functions that you expect to be common to
    many instruments should likely be coded there.
    '''

    telescope = 'APO'
    sitename = 'APO'

    # a string for the instrument name
    name = 'DIS'

    # are the slits in a "mask" (with different locations for every star)
    #               or a "longslit" (with one location for each star?)
    slitstyle = 'longslit'

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

    globallinekeys = ['AIRMASS', 'TELROT']

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
    keytosearch = 'FILENAME'

    # by what key should files be sorted in the summaries?
    summarysortkey = 'FILENAME'

    # within that header key, what words do we search for?
    wordstosearchfor = { 'dark':['dark'],
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


    def __init__(self, grating='R300'):
        '''
        This initializes a DIS instrument for a single camera side.
        '''

        # what grating is being used? ['B400', 'R300'] are supported so far
        self.grating = grating.upper()
        self.camera = self.grating[0]


        # you can use astropy to get the coordinates of an observatory
        # (we hard code them, so you can run this without internet,
        #  in case you might be trying to work on a plane or mountain top)
        # self.observatory = coord.EarthLocation.of_site(self.sitename)
        self.observatory = coord.EarthLocation.from_geodetic(-105.82*u.deg, 32.78*u.deg, 2798.0*u.m)


        # file search patten to get a *single* fileprefix for each exposure
        self.basicpattern = '*{}.fits'.format(self.camera.lower())

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
        self.namps = 1

        # FIXME -- it'd be real swell if these could be pulled from headers
        # a default, for now?
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


        # how many stitched images can we hold in memory?
        self.maximumimagesinmemory = 128
        self.maximumimagesinmemoryforscience = 32


    def setupDisperser(self):
        '''
        Setup the basics for the disperser.
        '''

        # what is the name and type of the disperser
        self.dispersertype = 'grating'
        self.disperser = self.grating

        # define a uniform grid of wavelengths for supersampling onto, later
        if self.grating == 'R300':
            self.uniformwavelengths = np.arange(4000, 10500)
            self.alignmentranges = {
                                                    r'$O_2$ - B':(6750,7050),
                                                    r'$O_2$ - A':(7500,7800),
                                                    r'Ca triplet':(8450,8750),
                                                    r'$H_2O$':(9200, 9700)
                                   }
            # pull good aligntment ranges from LDSS3C if you need them

        # the available arc lamps for wavelength calibration
        self.arclamps = ['He', 'Ne', 'Ar']

        # set up the wavelength calibration paths and files
        self.disperserDataDirectory = os.path.join(mosasaurusdirectory,
                                                'data/',
                                                self.name + '/',
                                                self.disperser + '/')


        self.wavelength2pixelsFile = os.path.join(self.disperserDataDirectory,
                '{0}_wavelength_identifications.txt'.format(self.grating))

        self.wavelengthsFile = os.path.join(self.disperserDataDirectory,
                'HeNeAr.txt')

        if self.binning == 2:
            self.offsetBetweenReferenceAndWavelengthIDs = 0# -1024
        else:
            self.offsetBetweenReferenceAndWavelengthIDs = 0
        # find the peak of the combined correlation function
        #if self.aperture.obs.instrument == 'LDSS3C':
        #    self.offsetBetweenReferenceAndWavelengthIDs = -1024 # KLUDGE KLUDGE KLUDGE! np.where(self.corre['combined'] == self.corre['combined'].max())[0][0] - len(x)
        #    # (old?) to convert: len(x) - xPeak = x + offsetBetweenReferenceAndWavelengthIDs
        #elif self.aperture.obs.instrument == 'IMACS': self.offsetBetweenReferenceAndWavelengthIDs = -75  # also a kludge



    def extractMidExposureTimes(self, headers):
        '''
        For an astropy table of extracted header values,
        extract the mid-exposure times in JD_UTC.

        Parameters
        ----------

        headers : astropy table (or generally iterable)
            a table containing temporal keywords, for determing times

        Returns
        -------

        times_earth : astropy time array
            the mid-exposure times, as measured at Earth

        '''

        # stitch together date+time strings
        timestrings = ['{0}'.format(row['DATE-OBS']) for row in headers]

        # calculate a JD from these times (and adding half the exposure time, assuming ut-time is the start of the exposure)
        starttimes = astropy.time.Time(timestrings, format='isot', scale='utc', location=self.observatory)

        # mid-exposure
        times_earth = starttimes + 0.5*headers['EXPTIME']*u.second

        # return as astropy times, in the UTC system, at the location of the telescope
        return times_earth


    def file2prefix(self, filename):
        '''
        This function returns a shortened fileprefix from a given filename.

        Parameters
        ----------

        filename : str
            The filename of a particular file.

        Returns
        -------

        prefix : str
            A shortened fileprefix for the file.
        '''
        tail = os.path.basename(filename)

        # let's pull out just the prefix from this DIS camera
        return tail.replace('{}.fits'.format(self.camera.lower()), '')

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
        return [prefix + '{}.fits'.format(self.camera.lower())]


    def gain(self, header):
        '''
        Return the gain, from a given header.
        (DIS has good headers, so this is easy.)
        '''
        return header['GAIN']

    def exptime(self, header):
        return header['EXPTIME']

    def darkexptime(self, header):
        return header['DARKTIME']

    def loadSingleCCD(self, filenames):
        '''
        Load a DIS image; subtract and trim its overscan.

        In general, this function should load and return
        a single image. If the detector uses multiple
        amplifiers to read out different parts of the
        same chip, this function should stitch those
        sections together.

        Parameters
        ----------
        filenames: list
            a list of relevant filenames (e.g. multiple amplifiers)

        Returns
        -------
        image: array
            a overscan-trimmed CCD image

        '''

        # for DIS, we need only one filename; it shouldn't be a list
        filename = filenames[0]

        # open (temporarily) the file
        with astropy.io.fits.open(filename) as hdu:
            self.data, self.header = hdu[0].data,  hdu[0].header

        # pull out just the data section
        dbottom, dtop, dleft, dright = iraf2python(self.header['DATASEC'])
        trimmed = self.data[dbottom:dtop, dleft:dright].astype(np.float)

         # subtract bias overscan (just one value)
        bbottom, btop, bleft, bright = iraf2python(self.header['BIASSEC'])
        biaslevel = np.median(self.data[bbottom:btop,bleft:bright])
        trimmed -= biaslevel

         # record that these have been stitched
        header = copy.copy(self.header)
        header['STITCHED'] = True

        # return the trimmed
        return trimmed.T, header
