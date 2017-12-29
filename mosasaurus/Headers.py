from .imports import *
import astropy.units as u, astropy.coordinates as coord
import astropy.table, astropy.time
from . import BJD

class Headers(Talker):

    '''An object to store the timeseries of image headers for this project -- good for keeping track of various external variables.'''
    def __init__(self, obs, **kwargs):
        '''Initialize a Headers object.'''

        # decide whether or not this creature is chatty
        Talker.__init__(self, **kwargs)

        # connect the observation object
        self.obs = obs

        # define a filename in which the header timeseries will be stored
        self.filename = self.obs.directory + 'headers.npy'

    def load(self, remake=True):
        '''Make sure the header table is loaded.'''

        self.speak('loading cube of image headers.')
        try:
            # does the headers attribute already exist?
            self.headers
            assert(remake == False)
            self.speak('header cube was already loaded')
        except (AssertionError, AttributeError):
            self.loadFromFile(remake=remake)


    def loadFromFile(self, remake=False):
        '''Load a table of header information from a pre-saved file.'''

        try:
            # try to load it from a pre-made file
            assert(remake==False)
            self.headers = astropy.table.Table(np.load(self.filename)[()])

            # make sure this table is the same length as the desired science exposures
            for k in self.headers.colnames:
                assert(len(self.headers['k']) == len(self.fileprefixes['science']))

            # say what happened
            self.speak('header timeseries loaded from {0}'.format(self.filename))

        except (AssertionError, IOError):
            # populate a header cube from the individual FITS headers
            self.loadFromScratch()

    def loadFromScratch(self):
       '''Extract parameters from science image headers, save as a structure.'''

       self.speak('looping through all science images to load their headers')

       # what keys do we need to populate?
       keys = self.obs.instrument.keysfortimeseries

       # create a dictionary of lists, to contain each key for all headers
       d = {}
       d['n'] = []
       for k in keys:
           d[k] = []

       # loop through the science images
       fileprefixes = self.obs.fileprefixes['science']
       for prefix in self.obs.fileprefixes:
           # get a number associated with this file
           n = self.obs.instrument.prefix2number(prefix)
           d['n'].append(n)

           # what is one file associated with this
           filename = self.obs.night.dataDirectory + self.instrument.prefix2files(prefix)[0]

           # load its header
           hdu = astropy.io.fits.open(filename)
           header = hdu[0].header
           for k in keys:
               d[k].append(header[k])

       # convert the dictionary of lists into a table
       self.headers = astropy.table.Table(d)[keys]
       print(self.headers)

       # convert times into more useful ones (including BJD)
       self.convertTimes()

       # save the table of headers
       np.save(self.filename, self.headers)
       self.speak('header cube saved to {0}'.format(self.filename))

    def convertTimes(self):
        '''Convert the header keys into more useful times; store them in the cube.'''
        self.speak('converting times from header into BJD')

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
