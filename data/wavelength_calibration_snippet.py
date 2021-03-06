import numpy as np
import matplotlib.pyplot as plt

# Here's how to plot your calibration lines
helium = np.load('/h/mulan0/data/working/K225b_ut161027_28/multipleapertures/aperture_605_1333/extracted_He.npy')[()]
neon = np.load('/h/mulan0/data/working/K225b_ut161027_28/multipleapertures/aperture_605_1333/extracted_Ne.npy')[()]
argon = np.load('/h/mulan0/data/working/K225b_ut161027_28/multipleapertures/aperture_605_1333/extracted_Ar.npy')[()]

plt.plot(helium['w'], helium[6.0]['raw_counts'])
plt.plot(neon['w'], neon[6.0]['raw_counts'])
plt.plot(argon['w'], argon[6.0]['raw_counts'])


# Here's how to read in the wavelength_identifications.txt file you made from matching up the lines
linelist = ascii.read('/h/mulan0/code/mosasaurus/data/150_wavelength_identifications.txt')
plt.plot(wavecalib['pixel'], wavecalib['wavelength'])
plt.plot(wavecalib['pixel'], wavecalib['wavelength'], 'bo')



# MISCELLANEOUS NOTES
# KEEP THE ARGON! IT'S LINED UP.
