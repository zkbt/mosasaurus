from imports import *
import astropy.units as u, astropy.coordinates as coord
import astropy.table, astropy.time
import BJD

class Headers(Talker):
    '''An object to store the timeseries of image headers for this project -- good for keeping track of various external variables.'''
    def __init__(self, obs, **kwargs):
        '''Initialize a Headers object.'''

        # decide whether or not this creature is chatty
        Talker.__init__(self, **kwargs)

        # add the observation object
        self.obs = obs
        self.filename = self.obs.workingDirectory + 'headers.npy'

    def load(self, remake=True):
        '''make sure the header table is loaded'''

        self.speak('loading cube of image headers')
        try:
            self.headers
            assert(remake == False)
            self.speak('header cube was already loaded')
        except:
            self.loadFromFile(remake=remake)


    def loadFromFile(self, remake=False):
        '''load a table of header information from a pre-saved file'''

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

       # what keys do we want to store?
       keys = ['date-obs', 'ut-date', 'ut-time', 'ut-end', 'scale', 'gain', 'epoch', 'airmass', 'ha', 'exptime', 'tempccd', 'templdss', 'focus', 'rotangle', 'rotatore']

       # create a dictionary of lists, to contain those for all headers
       d = {}
       d['n'] = []
       for k in keys:
           d[k] = []

       # loop through the science images
       ccdn = self.obs.nScience
       for n in self.obs.nScience:
           d['n'].append(n)
           filename = self.obs.dataDirectory+'ccd%04dc1.fits' % n
           hdu = astropy.io.fits.open(filename)
           header = hdu[0].header
           for k in keys:
               d[k].append(header[k])
           self.speak('   {0:10} {1:10} {2:10} {3:10}'.format(n, header['ut-date'],  header['ut-time'],  header['airmass']))

       # convert the dictionary of lists into a table
       self.headers = astropy.table.Table(d)

       # convert times into more useful ones (including BJD)
       self.convertTimes()

       # save the table of headers
       np.save(self.filename, self.headers)
       self.speak('header cube saved to {0}'.format(self.filename))

    def convertTimes(self):
        '''Convert the header keys into more useful times; store them in the cube.'''
        self.speak('converting times into BJD')


        # load one header, to get one-time information
        filename = self.obs.dataDirectory+'ccd%04dc1.fits' % self.obs.nScience[0]
        header = astropy.io.fits.open(filename)[0].header

        self.speak('loaded one header ({}) for site information'.format(filename))
        # stitch together date+time strings
        timestrings = ['{0} {1}'.format(row['ut-date'], row['ut-time']) for row in self.headers]

        # pull out the sitename
        sitename = header['SITENAME'].lower()
        observatory = coord.EarthLocation.of_site(sitename)

        # calculate a JD from these times (and adding half the exposure time, assuming ut-time is the start of the exposure)
        starttimes = astropy.time.Time(timestrings, format='iso', scale='utc', location=observatory)

        # mid-exposure
        times_earth = starttimes + 0.5*self.headers['exptime']*u.second

        # quote the JD_UTC
        self.headers['jd_utc']  = times_earth.jd

        # pull out RA and Dec (in degrees)
        ra, dec = header['RA-D'], header['DEC-D']

        # calculate BJD_TDB
        times_bary = BJD.toBJD(times_earth, ra=ra, dec=dec)
        self.headers['bjd'] = times_bary.jd
        self.headers['tdb-utc'] = times_earth.tdb.jd - times_earth.utc.jd
        self.headers['barycor'] = times_bary.tdb.jd - times_earth.tdb.jd
        
        a = self.input('Type enter if okay with BJD?!')
