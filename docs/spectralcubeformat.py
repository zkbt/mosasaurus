# this loads an example spectral cube (for WASP94Ab)
d = np.load('spectralCube_2stars_1626spectra.npy')[()]

# there are four different objects inside this dictionary
print d.keys()

# 1D arrays that correspond to the wavelength axis
print d['spectral']

# 1D arrays that correspond to the time axis
print d['temporal']

# 2D arrays for things that depend on star and time
#  with shape (nstars,nspectra)
print d['squares']

# 3D arrays for things that depend on star, time, and wavelength
#  with shape (nstars,nspectra,nwavelengths)
print d['cubes']
