from imports import *

class Night(Talker):
  '''Night objects handle information specific to the night on which an observation falls.'''
  def __init__(self, observation, **kwargs):
    '''Initialize a night object.'''
    self.obs = observation

    # decide whether or not this Reducer is chatty
    Talker.__init__(self, **kwargs)

  def obsLog(self, remake=False):
    '''Print and save an observing log, based on image headers.'''
    if self.obs.instrument == 'LDSS3C': keys = ['ut-time', 'object','aperture', 'grism', 'filter', 'exptime']
    elif self.obs.instrument == 'IMACS': keys = ['ut-time', 'object','slitmask', 'dispersr', 'filter', 'exptime']
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
