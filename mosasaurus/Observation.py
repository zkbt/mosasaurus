from .imports import *
from .Headers import Headers
#from Display import Display

#  an object that stores all the specifics related to a particular target/night of observing
class Observation(Talker):
    '''Observation object store basic information about an observation of
        one object with one instrument on one night.'''
    def __init__(self, target=None, instrument=None, night=None, **kwargs):
        '''Initialize an observation object.'''

        Talker.__init__(self)

        # set up connections to the other necessary objects
        self.target=target
        self.instrument=instrument
        self.night=night

        # make a directory hold all analyses for this observation
        self.directory = os.path.join(self.instrument.workingDirectory,
                                        "{}_{}".format(self.night.name, self.target.name))
        mkdir(self.directory)

        # set up the observation with the prefixes it will need
        try:
            self.setupFilePrefixes()
        except ValueError:
            self.speak('Hmmmmm...something funny happened with default file prefix choices. Please specify them by hand.')


    def __repr__(self):
        '''How should this object be represented as a string?'''
        return '[Observation of {} with {} on {}]'.format(self.target, self.instrument, self.night)

    def loadHeaders(self, remake=False):
        '''
        Load all the headers into one easy-to-manage table.
        '''

        h = Headers(self)
        h.load(remake=remake)
        self.headers = h.headers

    def setupFilePrefixes(self, **strategy):
        '''
        Define the image number arrays based on guesses from information
        in the file headers. These can be overwritten by custom setting the
        strategy attributes.

        For each file type (dark, bias, flat, science, reference, various arc lamps),
        you can specify one of ? ways in which to search for that kind of file

            None = fall back to the default, searching this instrument's
                    default header keyword for the default search string
        '''

        somethingisnew = False
        everything = (  self.instrument.detectorcalibrations +
                        self.instrument.arclamps +
                        self.instrument.extractables)

        self.exposures = {}
        choicesDirectory =  os.path.join(self.directory, 'files')
        # try
        for k in everything:
            fileofiles = os.path.join(choicesDirectory, 'filesfor{}.txt'.format(k))
            try:
                self.exposures[k] = astropy.io.ascii.read(fileofiles, delimiter='|')
                self.speak('loaded a list of [{}] files from {}'.format(k, fileofiles))
            except IOError:
                somethingisnew = True
                self.speak("couldn't find a list of files for [{}], making a guess".format(k))
                mkdir(choicesDirectory)

                # create a list of filenames
                if k in strategy.keys():

                    # if the strategy is a string, then just look for that
                    if type(strategy[k]) == list:
                        wordstosearchfor = strategy[k]

                    # add other options?

                else:
                    # find exposures where
                    if k == 'science':
                        wordstosearchfor = [self.target.starname]
                    elif k == 'reference':
                        wordstosearchfor = [self.target.starname]
                    else:
                        wordstosearchfor = self.instrument.wordstosearchfor[k]

                match = self.night.find(wordstosearchfor, self.instrument.keytosearch)
                self.exposures[k] = self.night.log[match]
                self.exposures[k].meta['comments'] = ['mosasaurus will treat these as [{}] exposures for {}'.format(k, self), '']
                self.exposures[k].write(fileofiles, **tablekw)
                self.speak('saved list of [{}] files to {}'.format(k, fileofiles))

        # give the user a chance to modify the initial (probably bad) guesses for which filenames to use
        if somethingisnew:
            answer = self.input("mosasaurus just created some guesses for filenames for {}".format(self) +
                                "\nPlease check the tables in {}".format(self.directory) +
                                "\nand edit them as necessary, before proceeding." +
                                "\n[Press enter to continue.]")

            self.speak('reloading from text files in {}, in case you made any changes'.format(choicesDirectory))
            for k in everything:
                fileofiles = os.path.join(choicesDirectory, 'filesfor{}.txt'.format(k))
                self.exposures[k] = astropy.io.ascii.read(fileofiles, delimiter='|')
                self.speak('loaded a list of [{}] files from {}'.format(k, fileofiles))

        self.fileprefixes = {}
        for k in everything:
            # pull out a simple array of filenames
            self.fileprefixes[k] = self.exposures[k]['fileprefix'].data

            # make the array indexable by the fileprefix
            self.exposures[k].add_index('fileprefix')


        # load a table of all the headers
        self.loadHeaders()
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
