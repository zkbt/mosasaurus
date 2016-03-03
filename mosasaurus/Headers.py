from imports import *
import astropy.table, astropy.time

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
        #for k in self.headers.keys():
        #    self.__dict__[k] = self.headers[k]

    def loadFromFile(self, remake=False):

        try:
            assert(remake==False)
            self.headers = astropy.table.Table(np.load(self.filename)[()])
            assert(len(self.headers['airmass']) == len(self.obs.nScience))
            self.speak('header cube loaded from {0}'.format(self.filename))
        except:
            self.loadFromScratch()

    def loadFromScratch(self):
       '''Extract parameters from science image headers, save as a structure.'''


       self.speak('looping through all science images to load their headers')
       keys = ['date-obs', 'ut-date', 'ut-time', 'ut-end', 'scale', 'gain', 'epoch', 'airmass', 'ha', 'exptime', 'tempccd', 'templdss', 'focus', 'rotangle', 'rotatore']
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

       self.headers = astropy.table.Table(d)
       self.convertTimes()
       np.save(self.filename, self.headers)
       self.speak('header cube saved to {0}'.format(self.filename))

    def convertTimes(self):
        '''Convert the header keys into more useful times; store them in the cube.'''
        self.speak('converting times into BJD (by calling IDL)')


        # stitch together date+time strings
        timestrings = ['{0} {1}'.format(row['ut-date'], row['ut-time']) for row in self.headers]

        # calculate a JD from these times (and adding half the exposure time, assuming ut-time is the start of the exposure)
        self.headers['jd']  = astropy.time.Time(timestrings, format='iso', scale='utc').jd + 0.5*self.headers['exptime']/24.0/60.0/60.0

        # use Jonathan Irwin's library to calculate BJD
        import BJD
        self.headers['bjd']  = [BJD.jdutc2bjdtdb(self.headers['jd'][i], self.obs.ra, self.obs.dec, observatory=self.obs.observatory) for i in range(len(self.headers['jd']))]
        print self.obs.ra, self.obs.dec
        self.headers['barycor'] = self.headers['bjd'] - self.headers['jd']
        print self.headers[['jd', 'barycor', 'bjd']]
        a = self.input('okay with BJD?!')
