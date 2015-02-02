from imports import *

class Headers(Talker):
    '''An object to store the timeseries of image headers for this project -- good for keeping track of various external variables.'''
    def __init__(self, obs, **kwargs):
        '''Initialize a Headers object.'''

        # decide whether or not this creature is chatty
        Talker.__init__(self, **kwargs)

        # add the observation object
        self.obs = obs
        self.filename = self.obs.workingDirectory + 'headers.npy'

    def load(self, remake=False):
        self.speak('loading cube of image headers')
        try:
            self.headers
            assert(remake == False)
            self.speak('header cube was already loaded')
        except:
            self.loadFromFile(remake=remake)
        for k in self.headers.keys():
            self.__dict__[k] = self.headers[k]

    def loadFromFile(self, remake=False):

        try:
            assert(remake==False)
            self.headers = np.load(self.filename)[()]
            assert(len(self.headers['airmass']) == len(self.obs.nScience))
            self.speak('header cube loaded from {0}'.format(self.filename))
        except:
            self.loadFromScratch()

    def loadFromScratch(self):
       '''Extract parameters from science image headers, save as a structure.'''


       self.speak('looping through all science images to load their headers')
       keys = ['date-obs', 'ut-time', 'ut-end', 'scale', 'gain', 'epoch', 'airmass', 'ha', 'exptime', 'tempccd', 'templdss', 'focus', 'rotangle', 'rotatore']
       d = {}
       d['n'] = []
       for k in keys:
           d[k] = []

       ccdn = self.obs.nScience
       for n in self.obs.nScience:
           d['n'].append(n)
           filename = self.obs.dataDirectory+'ccd%04dc1.fits' % n
           hdu = astropy.io.fits.open(filename)
           header = hdu[0].header
           for k in keys:
               d[k].append(header[k])
           self.speak('   {0:10} {1:10} {2:10} {3:10}'.format(n, header['ut-date'],  header['ut-time'],  header['airmass']))

       self.headers = d
       np.save(self.filename, self.headers)
       self.speak('header cube saved to {0}'.format(self.filename))

       '''ccdn = self.nScience
       date = [header['ut-date'] for header in headers]
       time = [header['ut-time'] for header in headers]
       airmass = [header['airmass'] for header in headers]
       timestrings = []
       for i in range(len(time)):
         timestrings.append( date[i] + ' ' + time[i])
       times = astropy.time.Time(timestrings, format='iso', scale='utc')
       np.savetxt(self.workingDirectory + 'headerInfo.txt', np.c_[ccdn, times.mjd, airmass])'''
