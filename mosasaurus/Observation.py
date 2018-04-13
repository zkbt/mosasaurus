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
        self.instrument.obs = self # link this observation back to the instrument
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

        # keep track of whether we need to re-load the fileprefixes
        somethingisnew = False

        # make a list of all the types of files we need to be aware of
        everything = (  self.instrument.detectorcalibrations +
                        self.instrument.arclamps +
                        self.instrument.extractables)

        # create a dictionary, whose keys will file types to keep track of
        self.exposures = {}

        # the directory where file choices should be saved
        choicesDirectory =  os.path.join(self.directory, 'files')


        # try loading existing files
        for k in everything:

            # the filename for fileprefixes for this type
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

                # this will find the indices that match the wordstosearchfor
                match = self.night.find(wordstosearchfor, self.instrument.keytosearch)

                # make sure this table is sorted by the fileprefix
                tableforthissubset = self.night.log[match]
                tableforthissubset.sort('fileprefix')


                self.exposures[k] = tableforthissubset
                self.exposures[k].meta['comments'] = ['mosasaurus will treat these as [{}] exposures for {}'.format(k, self), '']
                self.exposures[k].write(fileofiles, **tablekw)
                self.speak('saved list of [{}] files to {}'.format(k, fileofiles))

        # give the user a chance to modify the initial (probably bad) guesses for which filenames to use
        if somethingisnew:
            answer = self.input("mosasaurus just created some guesses for filenames for {}".format(self) +
                                "\nPlease check the tables in {}".format(choicesDirectory) +
                                "\nand edit them as necessary, before proceeding." +
                                "\nThis sets what exposures will be used for what."  +
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
