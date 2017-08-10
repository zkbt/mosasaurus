
from imports import *

class Night(Talker):
    '''Night objects handle information specific to the night.'''

    def __init__(self, name, instrument, **kwargs):
        '''Initialize a night object.'''

        Talker.__init__(self)

        # how do we refer to this night?
        self.name = name
        self.instrument = instrument

        # create an observing log for this night
        self.createNightlyLog()

    @property
    def dataDirectory(self):
        '''
        The absolute path to the directory containing the data for this night.
        '''
        return os.path.join(self.instrument.dataDirectory, self.name)


    def createNightlyLog(self):
        '''
        Print and save an observing log, based on image headers.
        '''

        # what are all the filenames associated with this night?
        pattern = '*.fits'

        self.filenames = glob.glob(os.path.join(self.dataDirectory, pattern))
        self.speak('Creating a log of all files on {}.'.format(self.name))
        self.speak('{} contains {} files matching "{}"'.format(
                            self.dataDirectory,
                            len(self.filenames),
                            pattern))

        # load or create a nightly obseravtion log
        self.logFilename = self.instrument.workingDirectory + 'nightly_log_{}.txt'.format(self.name)
        try:
            # load the nightly log as a astropy table
            self.log = astropy.io.ascii.read(self.logFilename)
            self.speak('Loaded a log file from {}.'.format(self.logFilename))
        except IOError:
            self.speak("Let's create a digital observing log from the image headers.")

            # print basic format of the first header
            self.speak("The format of the first file is:")
            hdu = astropy.io.fits.open(self.filenames[0])
            hdu.info()
            self.speak('')

            # extract the interesting keys from all the image headers
            rows = []
            for file in self.filenames:
                self.speak( 'Loading header information from {}.'.format(
                                os.path.basename(file)),
                            progress=True)

                keys, values = self.instrument.extractInterestingHeaderKeys(file)
                rows.append(dict(zip(keys, values)))

            # create an astropy table, and write it out to the working directory
            self.log = astropy.table.Table(data=rows, names=keys)
            
            self.log.write(self.logFilename,
                            format='ascii.fixed_width',
                            delimiter='|',
                            bookend=False)

    def find(self, wordstolookfor, placestolook):
        '''
        Find which rows of the log contain
        contain any one of the wordstolookfor
        in any one of the placestolook.

        (returns a boolean array)
        '''

        # create an array full of Falses
        match = np.zeros(len(self.log)).astype(np.bool)

        # find whether the wordstolookfor are seen anywhere in the placestolook
        for w in wordstolookfor:
            for p in placestolook:
                thismatch = np.array([w.lower() in word.lower() for word in self.log[p]])
                self.speak('{} elements in "{}" contained "{}"'.format(
                            np.sum(thismatch), p, w))
                match = match | thismatch

        # return the boolean array
        return match
