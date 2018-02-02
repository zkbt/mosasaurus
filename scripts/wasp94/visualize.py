from mosasaurus.Visualizer import Visualizer
from mosasaurus.imports import *
import glob


filenames = glob.glob('/Users/zkbt/Cosmos/Data/WASP94/*shifted*.npy')

f = filenames[1]
v = Visualizer(f)
if 'ut140801' in f:
    v.squishable.target='aperture_709_1066'
    v.squishable.comparisons=['aperture_751_1066']

if 'ut140805' in f:
    v.squishable.target='aperture_714_1064'
    v.squishable.comparisons=['aperture_755_1064']

if 'ut140809' in f:
    v.squishable.target='aperture_714_1063'
    v.squishable.comparisons=['aperture_755_1063']

#v.explore()


'''
z = v.squishable.corrected()
filtered = scipy.signal.medfilt(z, (3, 49))
cosmics = z - filtered
wavelengthstd = np.nanstd(cosmics, 0)
normalized = cosmics/wavelengthstd[np.newaxis, :]
notcosmics = np.abs(normalized) < 5
v.explore(normalized)
'''


wavelength = v.squishable.spectral['wavelength']
z = v.squishable.cubes['raw_counts'][v.squishable.target]
spectrum = np.nanmedian(z, 0)


remake=True
stride=50

'''
v.explore(key='raw_counts')
#v.overlayQuality()
v.loupe.moveCrosshair(y=np.min(wavelength), x=None)
v.loupe.ax['slicey'].set_xlabel('Photons/$\AA$')
v.loupe.ax['slicex'].set_ylabel('Photons/$\AA$')
plt.setp(v.loupe.ax['slicey'].get_xticklabels(), visible=False)
plt.savefig('raw.pdf', dpi=1000)
v.loupe.movieSlice(direction='y', filename='raw.mp4', remake=remake, stride=stride)
'''

v.explore(key='wavelengthed')
v.loupe.moveCrosshair(y=np.min(wavelength), x=None)
v.loupe.plotted['slicey'].set_data(spectrum, wavelength)
v.loupe.ax['slicey'].set_xlim(0, np.nanpercentile(spectrum, 99)*1.1)
plt.setp(v.loupe.ax['slicey'].get_xticklabels(), visible=False)
v.loupe.ax['slicex'].set_ylabel('Relative\nFlux')
plt.savefig('wavelengthed.pdf', dpi=1000)
v.loupe.movieSlice(direction='y', filename='wavelengthed.mp4', remake=remake, stride=stride)


'''
v.explore(key='corrected')
v.loupe.moveCrosshair(y=np.min(wavelength), x=None)
v.loupe.plotted['slicey'].set_data(spectrum, wavelength)
v.loupe.ax['slicey'].set_xlim(0, np.nanpercentile(spectrum, 99)*1.1)
plt.setp(v.loupe.ax['slicey'].get_xticklabels(), visible=False)
v.loupe.ax['slicex'].set_ylabel('Relative\nFlux')
plt.savefig('corrected.pdf', dpi=1000)
v.loupe.movieSlice(direction='y', filename='corrected.mp4', remake=remake, stride=stride)
'''


#v.explore(key='corrected')
#v.loupe.set_limits(0.98, 1.01)
#v.overlayQuality()
