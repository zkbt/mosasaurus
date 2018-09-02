# (see https://www.apo.nmsu.edu/arc35m/Instruments/DIS/images/henearr300w7700.gif)

# looking at some extracted spectra, this helps figure out the line list to use
# (I still have to add wavelengths by handing, but this saves some typing)
#cd /Users/zkbt/Cosmos/Data/APO/DIS/working/DIS-R300/UT180304/GJ649/extraction_default/aperture_390_0
import numpy as np
import craftroom.oned
from astropy.table import Table, vstack

elements = ['He', 'Ne', 'Ar']
tables = []
for e in elements:
    f = np.load('extracted_{}.npy'.format(e))[()]
    xPeak, yPeak, xfiltered, yfiltered = craftroom.oned.peaks(
                                                f['w'],
                                                f[8.0]['raw_counts'],
                                                plot=False,
                                                xsmooth=30,
                                                threshold=100,
                                                edgebuffer=10,
                                                widthguess=1,
                                                maskwidth=3,
                                                returnfiltered=True)
    N = len(xPeak)
    this = Table(data=[xPeak, yPeak, [' ']*N, [e]*N], names=['pixel', 'height', 'wavelength', 'name'])
    tables.append(this)
joined = vstack(tables)
joined.sort('pixel')
joined.write('/Users/zkbt/Dropbox/code/python/packages/mosasaurus/data/DIS/R300/R300-extracted.txt', delimiter='|', bookend=False, format='ascii.fixed_width')
