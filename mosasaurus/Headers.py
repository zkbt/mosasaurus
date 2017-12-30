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
        self.filename = os.path.join(self.obs.directory, 'headers.txt')

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
            #self.headers = astropy.table.Table(np.load(self.filename)[()])
            self.headers = astropy.io.ascii.read(self.filename, delimiter='|')

            # make sure this table is the same length as the desired science exposures
            for k in self.headers.colnames:
                assert(len(self.headers[k]) == len(self.obs.fileprefixes['science']))

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
       for prefix in fileprefixes:
           # get a number associated with this file
           n = self.obs.instrument.prefix2number(prefix)
           d['n'].append(n)

           # what is one file associated with this
           filename = os.path.join(self.obs.night.dataDirectory, self.obs.instrument.prefix2files(prefix)[0])

           # load its header
           hdu = astropy.io.fits.open(filename)
           header = hdu[0].header
           for k in keys:
               d[k].append(header[k])

       # convert the dictionary of lists into a table
       self.headers = astropy.table.Table(d)[keys]
       print(self.headers)

       # convert times into more useful ones (including BJD)
       self.determineBJD()

       # save the table of headers
       self.headers.write(self.filename, **tablekw)
       self.speak('header cube saved to {0}'.format(self.filename))

    def determineBJD(self):
        '''Convert the header keys into more useful times; store them in the cube.'''
        # (NOTE! This should eventually be moved to an instrument, because it'll vary from one to another)

        self.speak('converting times from header into BJD')

        # do the instrument-specific thing required to get a mid-transit time in JD_UTC
        times_earth = self.obs.instrument.extractMidExposureTimes(self.headers)

        # quote the JD_UTC in the header table
        self.headers['jd_utc']  = times_earth.jd

        # pull out RA and Dec (in degrees)
        icrs = self.obs.target.star.icrs

        # calculate BJD_TDB for one exposure, and print the results
        temp = BJD.toBJD(times_earth[0], icrs, verbose=True)

        # calculate BJD_TDB for all exposures
        times_bary = BJD.toBJD(times_earth, icrs, verbose=False)

        # and add them to the table (all in units of days)
        self.headers['bjd'] = times_bary.jd
        self.headers['tdb-utc'] = times_earth.tdb.jd - times_earth.utc.jd
        self.headers['barycor'] = times_bary.tdb.jd - times_earth.tdb.jd

        a = self.input("You just saw the calculation of the BJD times from information"
                        "\nin the FITS headers. It probably wouldn't be a bad idea to"
                        "\nspot-check this calculation against Jason Eastman's code at"
                        "\nhttp://astroutils.astronomy.ohio-state.edu/time/utc2bjd.html"
                        "\n[Press return if you're satisfied with the BJDs.]")
