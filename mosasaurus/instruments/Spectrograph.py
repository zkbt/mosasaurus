from ..imports import *

class Spectrograph(Talker):
    def __init__(self):

        Talker.__init__(self)
        self.setupDetector()
        self.setupDisperser()
        self.setupDirectories()
        self.setupExtraction()

    def extractInterestingHeaderKeys(self, file):
        '''
        This will work for files with only one interesting extension.
        If multiple extensions have interesting headers, this will
        need to be modified in the specific Spectrograph definition.
        '''

        # load the fits file
        hdu = astropy.io.fits.open(file)

        # for LDSS3C, one extension
        header = hdu[self.fitsextensionforheader].header

        # extract the values for these keys
        keys = self.keysforlogheader
        values = [header.get(k, '') for k in keys]

        # return a dictionary containing the interesting keys
        return keys, values
