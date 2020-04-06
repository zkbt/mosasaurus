from ..Spectrograph import *
from astropy.stats import mad_std
from astropy.coordinates import SkyCoord
import astropy.units as u
import illumination as il
from astropy.wcs import WCS
from astropy.io import fits


class WFC3(Spectrograph):
    '''
    This is an WFC3/IR instrument. It contains all instrument specific
    parameters and processes. Some procedures are inherited generally
    from Spectrograph, so functions that you expect to be common to
    many instruments should likely be coded there.
    '''

    telescope = 'HST'
    sitename = 'HST'

    # a string for the instrument name
    name = 'WFC3'

    # are the slits in a "mask" (with different locations for every star)
    #               or a "longslit" (with one location for each star?)
    slitstyle = 'slitless'

    # which header of the fits file contains the header with useful information
    fitsextensionforheader = 0

    # what keys should be included in the nightly logs for this instrument?
    # this is everything you might want to have access to at some point
    keysforlogheader = [    'ROOTNAME',
                            'FILENAME',
                            'INSTRUME',
                            'TARGNAME',
                            'RA_TARG',
                            'DEC_TARG',
                            'PROPOSID',
                            'SUNANGLE',
                            'MOONANGL',
                            'SUN_ALT',
                            'FGSLOCK',
                            'GYROMODE',
                            'DATE-OBS',
                            'TIME-OBS',
                            'EXPSTART',
                            'EXPEND',
                            'EXPTIME',
                            'PA_V3',
                            'POSTARG1',
                            'POSTARG2',
                            'OBSTYPE',
                            'SUBARRAY',
                            'SUBTYPE',
                            'FILTER',
                            'SAMP_SEQ',
                            'NSAMP',
                            'SAMPZERO',
                            'APERTURE',
                            'SAA_TIME',
                            'SCAN_TYP',
                            'SCAN_WID',
                            'ANG_SIDE',
                            'SCAN_ANG',
                            'SCAN_RAT',
                            'SCAN_LEN']

    # what keys should make it into condensed summary logs?
    # these are things you'll want for deciding what type each file is
    keysforsummary = [      'FILENAME',
                            'TARGNAME',
                            'EXPSTART',
                            'EXPTIME',
                            'FILTER',
                            'SUBTYPE',
                            'APERTURE']

    globallinekeys = [      'SUNANGLE',
                            'MOONANGL',
                            'SUN_ALT']

    # what keys do we want to store associated with a science timeseries?
    # these will show up, ultimately, in the 'temporal' key of a cube
    keysfortimeseries = [   'SUNANGLE',
                            'MOONANGL',
                            'SUN_ALT',
                            'FGSLOCK',
                            'DATE-OBS',
                            'TIME-OBS',
                            'EXPSTART',
                            'EXPEND',
                            'EXPTIME',
                            'PA_V3',
                            'SAMP_SEQ',
                            'NSAMP',
                            'SAMPZERO',
                            'SAA_TIME',
                            'SCAN_WID',
                            'ANG_SIDE',
                            'SCAN_ANG',
                            'SCAN_RAT',
                            'SCAN_LEN']


    # these keys are useful to search for guessing the filetype
    # for LDSS3C, the usful information is in the "object" key of the header
    # for other instruments, I could imagine "comments" being useful sometimes
    keytosearch = 'FILTER'



    # within that header key, what words do we search for?
    wordstosearchfor = { 'science':['G'],
                         'reference':['F']}

    wordstoavoid  =    { 'science':['F'],
                         'reference':['G']}
    # by what key should files be sorted in the summaries?
    summarysortkey = 'EXPSTART'

    def __init__(self, grism='G141'):
        '''
        This initializes a WFC3/IR instrument, with a particular grism.
        '''

        # what grism is being used (development stars with G141)
        self.grism = grism.upper()

        # let's pretend Hubble is at the center of the Earth, for simplicity
        self.observatory = coord.EarthLocation.from_geocentric(0*u.km, 0*u.km, 0*u.km)

        # file search patten to get a *single* fileprefix for each exposure
        self.basicpattern = '*ima.fits'

        # run the setup scripts, once these basics are defined
        Spectrograph.__init__(self)



    def setupDetector(self):
        '''
        Setup the basics for the detector.
        (These can be defaults can by modified by later changing
        the attributes, like for example:

            l = LDSS3C()
            l.gains = np.array([42.0, 54.0])
        '''

        # technically there are 4 amplifiers, but we'll pretend it's one
        self.namps = 1

        # WFC3 data is never binned
        self.binning = 1

        # what area of the detector contains real data?
        # (FIXME - this needs to be decided from the subarray)
        #self.dataleft = 0
        #self.dataright = 256
        #self.databottom = 0
        #self.datatop = 256

        # set up the size of the image
        self.xsize = np.inf#self.namps*(self.dataright - self.dataleft)
        self.ysize = np.inf#(self.datatop - self.databottom)

        # what are the calibrations we should expect
        self.detectorcalibrations = []


        # how many stitched images can we hold in memory?
        self.maximumimagesinmemory = 128
        self.maximumimagesinmemoryforscience = 32


    def setupDisperser(self):
        '''
        Setup the basics for the disperser.
        '''

        # what is the name and type of the disperser
        self.dispersertype = 'grism'
        self.disperser = self.grism

        # define a uniform grid of wavelengths for supersampling onto, later
        if self.disperser == 'G141':
            self.uniformwavelengths = np.arange(10500, 17000)
            self.alignmentranges = {
                                                    #r'$O_2$ - B':(6750,7050),
                                                    #r'$O_2$ - A':(7500,7800),
                                                    #r'Ca triplet':(8450,8750),
                                                    #r'$H_2O$':(9200, 9700)
                                                    # (maybe use the sodium line for M dwarfs?)
                                                    # and/or the water bandhead
                                   }
            # pull good aligntment ranges from LDSS3C if you need them

        # the available arc lamps for wavelength calibration
        self.arclamps = []

        # set up the wavelength calibration paths and files
        self.disperserDataDirectory = os.path.join(mosasaurusdirectory,
                                                    'data/',
                                                    self.name + '/',
                                                    self.disperser + '/')


        #self.wavelength2pixelsFile = os.path.join(self.disperserDataDirectory,
        #        '{0}_wavelength_identifications.txt'.format(self.grating))

        #self.wavelengthsFile = os.path.join(self.disperserDataDirectory,
        #        'HeNeAr.txt')

        #if self.binning == 2:
        #    self.offsetBetweenReferenceAndWavelengthIDs = 0# -1024
        #else:
        #    self.offsetBetweenReferenceAndWavelengthIDs = 0
        # find the peak of the combined correlation function
        #if self.aperture.obs.instrument == 'LDSS3C':
        #    self.offsetBetweenReferenceAndWavelengthIDs = -1024 # KLUDGE KLUDGE KLUDGE! np.where(self.corre['combined'] == self.corre['combined'].max())[0][0] - len(x)
        #    # (old?) to convert: len(x) - xPeak = x + offsetBetweenReferenceAndWavelengthIDs
        #elif self.aperture.obs.instrument == 'IMACS': self.offsetBetweenReferenceAndWavelengthIDs = -75  # also a kludge

#############
#START HERE!
#############


    def extractMidExposureTimes(self, headers):
        '''
        For an astropy table of extracted header values,
        extract the mid-exposure times in JD_UTC.

        Parameters
        ----------

        headers : astropy table (or generally iterable)
            a table containing temporal keywords, for determing times

        Returns
        -------

        times_earth : astropy time array
            the mid-exposure times, as measured at Earth

        '''

        # stitch together date+time strings
        timestrings = ['{0}'.format(row['DATE-OBS']) for row in headers]

        # calculate a JD from these times (and adding half the exposure time, assuming ut-time is the start of the exposure)
        starttimes = astropy.time.Time(timestrings, format='isot', scale='utc', location=self.observatory)

        # mid-exposure
        times_earth = starttimes + 0.5*headers['EXPTIME']*u.second

        # return as astropy times, in the UTC system, at the location of the telescope
        return times_earth


    def file2prefix(self, filename):
        '''
        This function returns a shortened fileprefix from a given filename.

        Parameters
        ----------

        filename : str
            The filename of a particular file.

        Returns
        -------

        prefix : str
            A shortened fileprefix for the file.
        '''
        tail = os.path.basename(filename)

        # let's pull out just the prefix from this DIS camera
        return tail.replace('_ima.fits', '')

    def prefix2number(self, prefix):
        '''
        This function returns a CCD number (not necessarily starting from 0)
        from a fileprefix.
        '''
        return 0 #np.int(prefix[-4:])

    def prefix2files(self, prefix):
        '''
        This function returns a list of filenames (without complete path)
        that are associated with this given prefix.
        '''
        return [prefix + '_ima.fits']


    def gain(self, header):
        '''
        Return the gain, from a given header.
        '''
        return 1.0

    def exptime(self, header):
        return header['EXPTIME']

    def darkexptime(self, header):
        return 0.0

    def loadSingleCCD(self, filenames):
        '''
        Load a WFC3 ima image.

        Subtract background from each read.
        (FIXME - This doesn't remove faint sources under ours!!)
        Stack the subtracted flux from individual reads together.
        The result should be a background-subtracted stacked image.

        Parameters
        ----------
        filenames: list
            a list of relevant filenames (e.g. multiple amplifiers)

        Returns
        -------
        image: array
            a overscan-trimmed CCD image

        '''


        filename = filenames[0]

        # FIXME -- move some of these to be instrument options!


        # open (temporarily) the file
        print(filename)
        with fits.open(filename) as hdu:


            # find all the SCI extensions inside the _ima
            n_extensions_per_read = 5
            nsamp = np.int((len(hdu)-1)/n_extensions_per_read)
            e_primary = 0
            e_sci = np.arange(nsamp)*n_extensions_per_read + 1
            print(e_sci)

            # transform the target location to pixel coordinates
            w = WCS(hdu[e_sci[0]].header)
            target_coord = SkyCoord(hdu[0].header['RA_TARG']*u.deg, hdu[0].header['DEC_TARG']*u.deg)
            xstart, ystart  = w.world_to_pixel(target_coord)
            print(xstart, ystart)

            # calculate the x and y pixel scales (in arcseconds/pixel)
            xscale = (np.sqrt(np.sum(w.pixel_scale_matrix[:,0]**2))*u.deg).to(u.arcsec)
            yscale = (np.sqrt(np.sum(w.pixel_scale_matrix[:,1]**2))*u.deg).to(u.arcsec)
            xscale, yscale

            print(xscale, yscale)
            # calculate the total scan height in pixels
            yscan_arcsec = hdu[0].header['SCAN_LEN']
            yscan_pixels = (yscan_arcsec*u.arcsec/yscale).decompose()

            # set up the basic scan geometry (FIXME -- make more flexible)
            yoffset = 0
            xoffset = 20

            boxwidth = 175
            boxleft = xstart + xoffset
            nboxes = hdu[0].header['NSAMP'] - 2
            buffer = 5
            yheight = yscan_pixels/(nboxes+1)
            boxheight = yheight*buffer



            # create some empty lists to populate with images
            rawreads = []
            maskedreads = []
            betterstacked = []
            stacked = []

            for i in range(nboxes):

                # pull out this science extension
                e = e_sci[i]

                # all WFC3 images have a 5 pixel border around them; trim that
                border = 5
                this_read = hdu[e].data[border:-border, border:-border]
                last_read = hdu[e+5].data[border:-border, border:-border]

                rows, cols = np.shape(this_read)
                x2d, y2d = np.meshgrid(np.arange(cols), np.arange(rows))

                boxcenter = ystart + yoffset + yheight*i #+ yheight/2.0
                boxbottom = boxcenter - boxheight/2.0
                boxright = boxleft + boxwidth
                boxtop = boxbottom + boxheight

                is_target = (y2d >= boxbottom)*(y2d <= boxtop)*(x2d >= boxleft)*(x2d <= (boxright))
                is_nottarget = is_target == False



                # pull out the flux accumuluated just in this read
                raw = this_read - last_read
                masked = raw*is_target

                #FIXME - still need to subtract background
                # estimate a crude background

                without_target = raw[is_nottarget]
                background_estimate = np.median(without_target)
                is_not_obvious_star = np.abs(without_target - background_estimate) < mad_std(without_target)*3
                better_background = np.median(without_target[is_not_obvious_star])

                raw_subtracted = raw - better_background
                masked = (raw - better_background)*is_target
                # estimate a better background and subtract it
                #
                #difference -= better_background

                #maskedread +
                rawreads.append(raw_subtracted)
                maskedreads.append(masked)
                betterstacked.append(np.sum(maskedreads,0))
                stacked.append(np.sum(rawreads, 0))


        movie_filename = filename.replace('fits', 'mp4')
        assert(movie_filename != filename)
        print(f'saving movie to {movie_filename}')


        frames = []

        # create a frame for the last - first image
        frames.append(il.imshowFrame(data=np.array(stacked), name='stacked'))

        # create a frame for each individual read
        frames.append(il.imshowFrame(data=np.array(rawreads), name='raw-reads'))

        # create a frame for each individual read, masked just around the target
        frames.append(il.imshowFrame(data=np.array(maskedreads),  name='masked-target'))

        # create a frame for each individual read, masked except for the target
        frames.append(il.imshowFrame(data=np.array(rawreads) - np.array(maskedreads), name='masked-other'))

        # create a frame for each individual read, masked except for the target
        frames.append(il.imshowFrame(data=np.array(betterstacked), name='better-stacked'))

        # make an illustration that contains all frames
        illust = il.GenericIllustration(imshows=frames, imshowrows=2)
        illust.plot()

        # pull out a single frame
        f = illust.frames['stacked']

        # plot the target reference location on the image
        plt.sca(f.ax)
        plt.scatter(xstart, ystart, c='black', marker='x', linewidth=0.5)
        plt.text(xstart, ystart, '\nTARG', ha='center', va='top', fontsize=6)


        # plot the individual-read boxes on each frame
        for i in range(nboxes):

            boxcenter = ystart + yoffset + yheight*i #+ yheight/2.0
            boxbottom = boxcenter - yheight/2.0

            for f in illust.frames.values():
                patch = plt.Rectangle([boxleft, boxbottom], boxwidth, yheight, fill=False)
                f.ax.add_patch(patch)

        illust.animate(filename=movie_filename, dpi=150)

        header = copy.copy(hdu[0].header)
        # return the trimmed
        return betterstacked[-1], header
