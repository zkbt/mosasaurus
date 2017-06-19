
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
            self.speak("Perhaps you didn't keep a perfect observing log; let's create one from image headers.")

            self.speak("The format of the first file is:")
            hdu = astropy.io.fits.open(self.filenames[0])
            hdu.info()
            self.speak('')

            rows = []
            for file in self.filenames:
                self.speak('Loading header information from {}.'.format(os.path.basename(file)), progress=True)

                keys, values = self.instrument.extractInterestingHeaderKeys(file)
                rows.append(dict(zip(keys, values)))

            self.log = astropy.table.Table(data=rows, names=keys)
            self.log.write(self.logFilename,
                            format='ascii.fixed_width',
                            delimiter='|',
                            bookend=False)
'''
  def obsLog(self, remake=False):


    logFile = self.obs.workingDirectory + 'complete_observing_log.txt'
    if os.path.exists(logFile) and remake==False:
      self.speak('looks like an observing log exists at {0}'.format(logFile))
    else:
      self.speak("perhaps you didn't keep a perfect observing log; let's create one from image headers")
      log = []
      fileprefixes = [x.split('c1')[0] for x in glob.glob(self.obs.dataDirectory + '*c1.fits')]
      for fileprefix in fileprefixes:
        hdu = astropy.io.fits.open(fileprefix+'c1.fits')

        string = "{0: <9}".format(fileprefix.split('/')[-1])
        for k in keys:
          string += ' ' + zachopy.utils.truncate(str(hdu[0].header[k]),n=12 + 3*(k == 'object'))
        self.speak(string)
        log.append(string)
      f = open(logFile, 'w')
      f.writelines(log)
      f.close
      self.speak("a digital observing log saved to {0}".format(logFile))
'''
