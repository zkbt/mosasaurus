
from .imports import *

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


    def __repr__(self):
        '''How should this object be represented as a string?'''
        return '<Night {}>'.format(self.name)

    @property
    def dataDirectory(self):
        '''
        The absolute path to the directory containing the data for this night.
        '''
        return os.path.join(self.instrument.dataDirectory, self.name)


    def createNightlyLog(self, remake=False):
        '''
        Print and save an observing log, based on image headers.
        '''

        # what are all the filenames associated with this night?
        pattern = self.instrument.basicpattern

        self.filenames = glob.glob(os.path.join(self.dataDirectory, pattern))
        self.speak('Creating a log of all files on {}.'.format(self.name))
        self.speak('{} contains {} files matching "{}"'.format(
                            self.dataDirectory,
                            len(self.filenames),
                            pattern))

        # load or create a nightly obseravtion log
        self.logFilename = self.instrument.workingDirectory + 'nightly_log_{}.txt'.format(self.name)

        try:
            assert(remake == False)
            # load the nightly log as a astropy table
            self.log = astropy.io.ascii.read(self.logFilename)
            self.speak('Loaded a log file from {}.'.format(self.logFilename))
        except (AssertionError, IOError):
            self.speak("Let's create a digital observing log from the image headers.")

            # print basic format of the first header
            self.speak("The format of the first file is:")
            hdu = astropy.io.fits.open(self.filenames[0])
            hdu.info()
            self.speak('')

            # extract the interesting keys from all the image headers
            self.rows = []
            for file in self.filenames:
                self.speak( 'Loading header information from {}.'.format(
                                os.path.basename(file)), progress=True)

                try:
                    keys, values = self.instrument.extractInterestingHeaderKeys(file)
                except:
                    self.speak('There was something troubling about {}'.format(file))
                    continue

                for i, v in enumerate(values):
                    if type(v) == astropy.io.fits.card.Undefined:
                        values[i] = '???'
                this = dict(zip(keys, values))

                # record the fileprefix associated with this
                this['fileprefix'] = self.instrument.file2prefix(file)
                # clean up comments with lots of newlines (mysteriously a problem at least for LDSS3C)
                try:
                    this['comment'] = str(this['comment']).strip()
                except KeyError:
                    pass


                self.rows.append(this)

            # create an astropy table, and write it out to the working directory
            self.log = astropy.table.Table(data=self.rows, names=['fileprefix'] + keys)

            self.log.write(self.logFilename,
                            format='ascii.fixed_width',
                            delimiter='|',
                            bookend=False, overwrite=remake)
            self.speak('wrote a log for {} to {}'.format(self.name, self.logFilename))

        # create a condensed summary of the log
        self.createSummaryLog(remake=remake)

    def createSummaryLog(self, remake=False):
        '''
        This function helps someone running the code figure out
        what unique types of exposures exist, so they can pick
        which ones to use for calibrations, for references, and
        for science.
        '''

        # try to not to duplicate effort
        self.summaryFilename = self.instrument.workingDirectory + 'nightly_summary_{}.txt'.format(self.name)
        try:
            assert(remake == False)
            # load the nightly summary as a astropy table
            self.summarylog = astropy.io.ascii.read(self.summaryFilename)
            self.speak('Loaded a condensed summary from {}.'.format(self.summaryFilename))
        except (AssertionError, IOError):
            self.speak("Let's summarize the nightly log.")

            # group the nightly log by the sorting key
            column = self.instrument.keytosearch


            # fix any that have ill-defined entries
            if self.log.mask != None:
                bad = (self.log[column].mask == True)
                self.log[column][bad] = '???'
                # (this should break if the column's not a string)
            bad = (self.log[column] == '')|np.array([type(c) == None for c in self.log[column] ])
            self.log[column][bad] = '???'

            # group by that sorting column
            subset = self.log[self.instrument.keysforsummary]
            grouped = (subset.group_by(column)).copy()


            # aggregrate the rows together
            warnings.filterwarnings('ignore')
            def glom(a):
                try:
                    return '{:.3}\pm{:.3}'.format(np.nanmedian(a), 1.486*np.median(np.abs(a - np.median(a))))
                except TypeError:
                    s = ' + '.join(['("{}")x{}'.format(u, np.sum(a == u)) for u in np.unique(a)])
                    if len(s) > 30:
                        return '("{}") etc...'.format(a[0])
                    else:
                        return s

            # create and write an aggregrated table
            self.summarylog = grouped.groups.aggregate(glom)
            count = [np.sum(self.log[column] == c) for c in self.summarylog[column]]
            self.summarylog['count'] = count
            self.summarylog = self.summarylog[['count'] + grouped.colnames]

            # sort this by the file prefixes, again
            self.summarylog.sort(self.instrument.summarysortkey)


            self.summarylog.write(self.summaryFilename, **tablekw)

            self.speak('wrote the aggregated summary for {} to {}'.format(self.name, self.summaryFilename))


    def find(self, wordstolookfor, placetolook, wordstoavoid=[]):
        '''
        Find which rows of the log contain
        contain any one of the wordstolookfor
        in any one of the placetolook.

        Parameters
        ----------

        wordstolookfor : list of strings
            If a string matches any of these strings,
            it will be included in the list of returned rows.

        placestolook : string
            Which column of the log should be searched?

        wordstoavoid : list of strings
            If a string matches any of these strings,
            it will *not* be included in the list of returned rows.

        (returns a boolean array)
        '''

        # create an array full of Falses
        match = np.zeros(len(self.log)).astype(np.bool)

        # find whether the wordstolookfor are seen anywhere in the placetolook
        for w in wordstolookfor:
            thismatch = np.array([w.lower() in word.lower() for word in self.log[placetolook]])

            # reject those bad ones
            for bad in wordstoavoid:
                isbad = np.array([bad.lower() in word.lower() for word in self.log[placetolook]])
                thismatch = thismatch & (isbad == False)

            self.speak('{} elements in "{}" contained "{}"'.format(
                        np.sum(thismatch), placetolook, w))
            match = match | thismatch

        # return the boolean array
        return match
