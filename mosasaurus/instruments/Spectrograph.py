from ..imports import *

def iraf2python(s):
    '''
    Convert an IRAF (1-indexed) column-row string ('[c1:c2,r1:r2]')
            to Python (0-indexed) [r1:r2,c1:c2]

    Returns
    -------
    s : str
        an IRAF-style column-row string like '[c1:c2,r1:r2]''

    bottom, top, left, right : int
        The corners of an array to extract, with Python indices.

    '''
    cols, rows = s.strip('[]').split(',')
    bottom, top = np.array(rows.split(':')).astype(np.int) - np.array([1,0])
    left, right = np.array(cols.split(':')).astype(np.int) - np.array([1,0])
    return bottom, top, left, right


class Spectrograph(Talker):

    # what's the basic path where *everything* will be stored?
    path = '/Users/zkbt/Cosmos/Data/'
    # FIXME -- should this *try* to pull from an environment variable?


    # are the slits in a "mask" (with different locations for every star)
    #               or a "longslit" (with one location for each star?)
    slitstyle = 'mask'

    def __init__(self, **kw):
        '''
        Initialize a spectrograph, using custom-defined
        Detector, Disperser, Directories, and Extractions
        that are specified for the individual instruments.

        Parameters
        ----------

        **kwargs will be passed to *every* setup function
        '''

        Talker.__init__(self)
        self.setupDetector(**kw)
        self.setupDisperser(**kw)
        self.setupDirectories(**kw)
        self.setupExtraction(**kw)

    def __repr__(self):
        '''
        How should this object be represented as a string?
        '''
        return '<Spectrograph {}-{}>'.format(self.name, self.disperser)

    def extractInterestingHeaderKeys(self, file):
        '''
        This will work for files with only one interesting extension.
        If multiple extensions have interesting headers, this will
        need to be modified in the specific Spectrograph definition.

        Parameters
        ----------
        file : string
            Filename to open, and extract some keys from the header.
        '''

        # load the fits file
        hdu = astropy.io.fits.open(file)

        # for LDSS3C, one extension
        header = hdu[self.fitsextensionforheader].header

        # extract the values for these keys
        keys = self.keysforlogheader
        values = [header.get(k, None) for k in keys]

        # return a dictionary containing the interesting keys
        return keys, values

    @property
    def zapcosmics(self):
        return self.extractiondefaults['zapcosmics']

    def setupExtraction(self):
        '''
        Setup the default extraction parameters associated with this instrument.

        (FIXME? This sets up defaults, which are then actually used. Perhaps
         include a step where these are actually adopted? Or, just replace
         everywhere with "extractionparameters"?
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

    @property
    def directory():
        '''Kludgy shortcut for this instrument's working directory. FIXME?'''
        return self.workingDirectory

    def setupDirectories(self,
            baseDirectory=None,
            dataDirectory='data/',
            workingDirectory='working/',
            extractionDirectory='extraction/'):
        '''
        Setup the basic file structure associated with this instrument.
        If you store your data anywhere other than the defaults listed here,
        then adjust the directory structure with:
            i = LDSS3C()
            i.setupDirectories("/absolute/path/to/base/directory/")

            $MOSASAURUSDATA/

        '''

        # absolute path to where the data are located
        if baseDirectory is None:
            # default to the telescope name + instrument name as directories
            baseDirectory = os.path.join(self.telescope, self.name)
        self.baseDirectory = os.path.join(self.path, baseDirectory)
        self.speak('set the base directory to {}'.format(self.baseDirectory))

        # data directory is where raw (not-to-be-modified) files are stored
        self.dataDirectory = os.path.join(self.baseDirectory,
                                                dataDirectory)
        mkdir(self.dataDirectory)
        self.speak('expecting to find raw data in nightly directories inside {}'.format(self.dataDirectory))



        # working directory is where in-progress files + results are created
        self.spectrographWorkingDirectory = os.path.join(self.baseDirectory,
                                                  workingDirectory)
        mkdir(self.spectrographWorkingDirectory)
        self.speak('all working data will be stored in {}'.format(self.spectrographWorkingDirectory))

        # extraction directory is where extracted spectra land
        self.spectrographExtractionDirectory = os.path.join(self.baseDirectory,
                                                            extractionDirectory)
        mkdir(self.spectrographExtractionDirectory)

        # create directories
        # directory for data associated generally with the disperser
        self.workingDirectory = os.path.join(self.spectrographWorkingDirectory,
                                             '{}-{}'.format(self.name, self.disperser))
        mkdir(self.workingDirectory)

        self.extractionDirectory = os.path.join(self.spectrographExtractionDirectory,
                                                self.disperser)
        mkdir(self.extractionDirectory)
