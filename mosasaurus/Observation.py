from .imports import *
from .Headers import Headers
#from Display import Display

#  an object that stores all the specifics related to a particular target/night of observing
class Observation(Talker):
    '''Observation object store basic information about an observation of one object on one night.'''
    def __init__(self, target=None, instrument=None, night=None, **kwargs):
        '''Initialize an observation object.'''

        Talker.__init__(self)

        # set up connections to the other necessary objects
        self.target=target
        self.instrument=instrument
        self.night=night

        # come up with a guess for the filenames
        # make a guess, check with user they're good, then write a file saying they're confirmed

        #self.loadHeaders()

        #self.fileprefixes = self.fileprefix(self.nNeeded)

    def __repr__(self):
        '''How should this object be represented as a string?'''
        return '[Observation of {} with {} on {}]'.format(self.target, self.instrument, self.night)

    def loadHeaders(self, remake=False):
        self.headers = Headers(self, mute=self._mute, pithy=self._pithy)
        self.headers.load(remake=remake)

    def guessFiles(self):
        '''
        Define the image number arrays based on guesses from information
        in the file headers. These can be overwritten by custom setting the
        self.n* attributes.
        '''

        self.instrument.findBiases(self.night)

        '''
        # modify these to make some guesses -- will need to know the mask name for science data
        self.nReference = np.arange(int(dictionary['nReference'][0]), int(dictionary['nReference'][1])+1)
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
        '''



        self.cosmicThreshold = float(dictionary['cosmicThreshold'])
        self.nNeeded = np.concatenate((self.nReference, self.nScience, self.nHe, self.nNe, self.nAr, self.nWideFlat, self.nWideMask, self.nThinMask, self.nFinder))
        self.cal_dictionary = {'He':self.nHe, 'Ne':self.nNe, 'Ar':self.nAr, 'Reference':self.nReference, 'WideFlat':self.nWideFlat, 'WideMask':self.nWideMask, 'ThinMask':self.nThinMask, 'Bias':self.nBias, 'Dark':self.nDark, 'Science':self.nScience, 'Finder':self.nFinder}




        self.ra =np.float(dictionary['ra'])
        self.dec =np.float(dictionary['dec'])
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
