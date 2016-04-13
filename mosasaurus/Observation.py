from imports import *
from Headers import Headers
#from Display import Display

#  an object that stores all the specifics related to a particular target/night of observing
class Observation(Talker):
    '''Observation object store basic information about an observation of one object on one night.'''
    def __init__(self, filename, nods9=False, **kwargs):
        '''Initialize an observation object.'''

        # decide whether or not this creature is chatty
        Talker.__init__(self, **kwargs)


        self.readParameters(filename)

        self.fileprefixes = self.fileprefix(self.nNeeded)

        zachopy.utils.mkdir(self.workingDirectory)
        #self.display = Display(nods9=nods9)


    def loadHeaders(self, remake=False):
        self.headers = Headers(self, mute=self.mute, pithy=self.pithy)
        self.headers.load(remake=remake)

    def readParameters(self, filename):
        '''A function to read in a stored parameter file, with all details needed for extraction.'''
        self.speak('trying to read {0} for observation parameters'.format(filename))
        file = open(filename)
        lines = file.readlines()
        dictionary = {}
        for i in range(len(lines)):
          if lines[i] != '\n' and lines[i][0] != '#':
            split = lines[i].split()
            key = split[0]
            entries = split[1:]
            if len(entries) == 1:
              entries = entries[0]
            dictionary[key] = entries
        self.name = dictionary['name']
        self.night = dictionary['night']
        self.grism = dictionary['grism'].lower()
        self.instrument = dictionary['instrument']
        if "LDSS" in self.instrument:
            self.observatory = 'lco'
        self.baseDirectory = dictionary['baseDirectory']

        # set up the wavelength calibration paths
        self.referenceDirectory = mosasaurusdirectory + 'data/'
        self.wavelength2pixelsFile = self.referenceDirectory  + '{0}_wavelength_identifications.txt'.format(self.grism)
        self.wavelengthsFile = self.referenceDirectory + 'HeNeAr.txt'

        zachopy.utils.mkdir(self.baseDirectory + dictionary['workingDirectory'])
        self.workingDirectory = self.baseDirectory + dictionary['workingDirectory'] + self.name + '_' + self.night +'/'
        zachopy.utils.mkdir(self.workingDirectory)
        self.dataDirectory = self.baseDirectory + dictionary['dataDirectory'] + self.night +'/'
        zachopy.utils.mkdir(self.dataDirectory)
        self.extractionDirectory = self.workingDirectory + dictionary['extractionDirectory']
        zachopy.utils.mkdir(self.extractionDirectory)
        self.extractionWidth = int(dictionary['extractionWidth'])
        self.skyGap = int(dictionary['skyGap']    )
        self.skyWidth = int(dictionary['skyWidth'])
        self.cosmicThreshold = float(dictionary['cosmicThreshold'])
        self.nUndispersed = np.arange(int(dictionary['nUndispersed'][0]), int(dictionary['nUndispersed'][1])+1)
        self.nScience = np.arange( int(dictionary['nScience'][0]),  int(dictionary['nScience'][1])+1)
        self.nHe = np.arange( int(dictionary['nHe'][0]),  int(dictionary['nHe'][1])+1)
        self.nNe = np.arange( int(dictionary['nNe'][0]),  int(dictionary['nNe'][1])+1)
        self.nAr = np.arange( int(dictionary['nAr'][0]),  int(dictionary['nAr'][1])+1)
        self.nDark = np.arange( int(dictionary['nDark'][0]),  int(dictionary['nDark'][1])+1)
        self.nWideFlat = np.arange(int(dictionary['nWideFlat'][0]),  int(dictionary['nWideFlat'][1])+1)
        self.nWideMask = np.arange(int(dictionary['nWideMask'][0]),  int(dictionary['nWideMask'][1])+1)
        self.nThinMask = np.arange(int(dictionary['nThinMask'][0]),  int(dictionary['nThinMask'][1])+1)
        if len(dictionary['nFinder']) == 1:
          self.nFinder = np.array([int(dictionary['nFinder'])])
        else:
          self.nFinder = np.arange(int(dictionary['nFinder'][0]),  int(dictionary['nFinder'][1])+1)
        self.nBias = np.arange(int(dictionary['nBias'][0]),  int(dictionary['nBias'][1])+1)
        self.nNeeded = np.concatenate((self.nUndispersed, self.nScience, self.nHe, self.nNe, self.nAr, self.nWideFlat, self.nWideMask, self.nThinMask, self.nFinder))
        self.cal_dictionary = {'He':self.nHe, 'Ne':self.nNe, 'Ar':self.nAr, 'Undispersed':self.nUndispersed, 'WideFlat':self.nWideFlat, 'WideMask':self.nWideMask, 'ThinMask':self.nThinMask, 'Bias':self.nBias, 'Dark':self.nDark, 'Science':self.nScience, 'Finder':self.nFinder}
        self.traceOrder    =  int(dictionary['traceOrder'])
        self.widthGuess = int(dictionary['widthGuess'])
        self.nFWHM    = float(dictionary['nFWHM'])
        self.blueward = int(dictionary['blueward'])
        self.redward = int(dictionary['redward'])
        self.target = [int(x) for x in dictionary['target']]
        self.goodComps = [int(x) for x in dictionary['comp']]
        self.dataleft    = int(dictionary['dataleft'])
        self.dataright    = int(dictionary['dataright'])
        self.databottom    = int(dictionary['databottom'])
        self.datatop    = int(dictionary['datatop'])
        self.namps =int(dictionary['namps'])
        self.xsize = self.namps*(self.dataright - self.dataleft)
        self.ysize = (self.datatop - self.databottom)
        self.ra =np.float(dictionary['ra'])
        self.dec =np.float(dictionary['dec'])
        self.binning = np.float(dictionary['binning'])
        self.subarray = np.float(dictionary['subarray'])
        try:
            self.gains = [float(x) for x in dictionary['gains']]
        except (ValueError,KeyError):
            self.gains = None

        try:
            self.slow = bool(dictionary['slow'])
        except KeyError:
            self.slow = False

        self.correlationAnchors = [float(x) for x in dictionary['correlationAnchors']]
        self.correlationRange = [float(x) for x in dictionary['correlationRange']]
        self.correlationSmooth = float(dictionary['correlationSmooth'])
        self.cosmicAbandon = float(dictionary['cosmicAbandon'])
        self.speak('observation parameters have been read and stored'.format(filename))

        self.displayscale=0.25
    def fileprefix(self, n):
        '''Feed in a ccd number, spit out the file prefix for that CCD amplifier pair.'''
        try:
          return [self.dataDirectory + 'ccd{0:04}'.format(x) for x in n]
        except:
          return self.dataDirectory + 'ccd{0:04}'.format(n)
