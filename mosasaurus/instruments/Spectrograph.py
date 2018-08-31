from ..imports import *

def iraf2python(s):
    '''
    Convert an IRAF (1-indexed) column-row string ('[c1:c2,r1:r2]')
            to Python (0-indexed) [r1:r2,c1:c2]

    Returns
    -------
    s : str
        an IRAF-style column-row string like '[c1:c2,r1:r2]''

    bottom, top, left, right : int
        The corners of an array to extract, with Python indices.

    '''
    cols, rows = s.strip('[]').split(',')
    bottom, top = np.array(rows.split(':')).astype(np.int) - np.array([1,0])
    left, right = np.array(cols.split(':')).astype(np.int) - np.array([1,0])
    return bottom, top, left, right


class Spectrograph(Talker):
    def __init__(self):

        Talker.__init__(self)
        self.setupDetector()
        self.setupDisperser()
        self.setupDirectories()
        self.setupExtraction()

    def __repr__(self):
        '''
        How should this object be represented as a string?
        '''
        return '<Spectrograph {}>'.format(self.name)

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
        values = [header.get(k, None) for k in keys]

        # return a dictionary containing the interesting keys
        return keys, values

    @property
    def zapcosmics(self):
        return self.extractiondefaults['zapcosmics']
