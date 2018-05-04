# Elisabeth and Hannah work things out

# from GJ 1132/working/ 
# wavecal[0] is an N-degree Legenre polynomial; wavecal[1] is the range over which it works
wavecal = np.load('aperture_832_1173_wavelengthcalibration.npy')[()]

# from GJ 1132/working/
stamp = np.load('calibStamps_aperture_832_1173.npy')[()

# recreate the conversion from pixel position to wavelength
from numpy.polynomial.legendre import Legendre 
L = Legendre(wavecal[0], domain=wavecal[1])

# make an interpolation to go from wavelength to pixel position
# stamp['w'] is the pixel range centered on the star, e.g., array([-1173, -1172, -1171, ...,   872,   873,   874])

def getPixelPositon (wavelength):                                   
    return np.interp(wavelength, L(stamp['w'][:,0]), stamp['w'][:,0])   




# some wavelengths from LDSS3C manual calibration; WE NEED TO MAKE THESE FOR IMACES GRISMS!
from astropy.io import ascii
waves = ascii.read('/h/mulan0/code/mosasaurus/data/IMACS/gri-150-10.8/gri-150-10.8_wavelength_identifications.txt')
