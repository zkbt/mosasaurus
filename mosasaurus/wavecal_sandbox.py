import chromatic.starspots.phoenix as p
from scipy.signal import medfilt
import numpy as np
import matplotlib.pyplot as plt
import rainbowconnection
from craftroom.resample import bintogrid
import astropy.units as u

import tqdm

def find_wavecal(e, aperture=8.0, Teff=7000, N = 50, downsample = 3, 
                 w0=[4700, 5000], dwdp=[1.0, 5.0], filtpixels=100):
    
    # read the spectrum
    wnm, sinspace = p.read_phoenix(Teff=Teff, R=2000)
    w = wnm*10 # convert to AA
    
    earth = rainbowconnection.Earth()
    s = sinspace*earth.transmission(wnm*u.nm)
    
    plt.plot(w, sinspace)
    plt.plot(w, s)
    plt.xlim(3000, 8000)
    plt.show()

    uw, uf = e['w'], e[aperture]['raw_counts']
    totalphotons = np.maximum(e[aperture]['raw_counts'] + e[aperture]['sky']*16, 100)
    ue = np.sqrt(totalphotons)

    dw, df, de = bintogrid(uw, uf, ue, newx=uw[::downsample])

    medscale = int(np.round(filtpixels/downsample/2)*2 + 1)
    dsmoothed = medfilt(df, medscale)


    def model(a=5000, b=5):
        wave = a + b*dw #+ c*(pix - opix)**2
        nw, nf = bintogrid(w, s, newx=wave, drop_nans=False)
        nsmoothed = medfilt(nf, medscale)
        #plt.plot(dw, dsmoothed)
        blaze_approx = dsmoothed/nsmoothed
        blaze_approx[np.isfinite(blaze_approx) == False] = 0

        #plt.plot(dw, df, alpha=0.5)
        #plt.plot(dw, nf*blaze_approx, alpha=0.5)
        #plt.plot(dw, nsmoothed)
        #plt.plot(dw, blaze_approx*1e4)
        return nw, nf*blaze_approx
    #plt.plot(dw, df)

    def residuals(a=5000, b=5, visualize=False):

        nw, nf = model(a=a, b=b)
        nf[np.isfinite(nf) == False] = 0.0
        r = (df - nf)/de
        chisq = np.sum(r**2)
        if visualize:
            plt.plot(dw, df, alpha=0.5)
            plt.plot(dw, nf, alpha=0.5)
            plt.title(f'a={a:.4}, b={b:.4}, chisq={chisq:.3}')

        return chisq


    # set up a grid
    a1d = np.linspace(w0[0], w0[1], N)
    b1d = np.linspace(dwdp[0], dwdp[1], N)

    # do a grid search
    chisq = np.zeros((N, N))
    for i, a in tqdm.tqdm(enumerate(a1d)):
        for j, b in enumerate(b1d):
            chisq[i,j] = residuals(a=a, b=b)

    # pick out the best
    ibest, jbest = np.where(chisq == np.min(chisq))
    abest = a1d[ibest[0]]
    bbest = b1d[jbest[0]]

    fi, ax = plt.subplots(1, 2)
    
    # plot the results
    plt.sca(ax[1])
    plt.imshow(np.log10(chisq),
               origin='lower',
               extent=[np.min(b1d), np.max(b1d), np.min(a1d), np.max(a1d)],
              aspect='auto'); 
    plt.colorbar(label='log($\chi^2$)')
    plt.scatter(bbest, abest, edgecolor='white', facecolor='none', marker='o', lw=3, s=500)
    plt.xlabel('b = d$\lambda$/dw')
    plt.ylabel('a = $\lambda$ at w=0');

    plt.sca(ax[0])
    residuals(a=abest, b=bbest, visualize=True)
    return  bbest, abest


from mosasaurus.WavelengthCalibrator import WavelengthCalibrator, Talker
class PHOENIX_WavelengthCalibrator(WavelengthCalibrator):
    def __init__(self,  aperture,
                        polynomialdegree=1, **kw):
        Talker.__init__(self)

        # keep track of the aperture this belongs to
        self.aperture = aperture

        # set up the some defaults
        self.polynomialdegree = polynomialdegree

        # either load a previous wavecal, or create a new one
        self.populate(**kw)
        
    def populate(self, restart=False, **kw):
        '''
        Make sure the wavelength calibration is populated,
        either by loading it or making it.
        '''
        try:
            # populate the wavelength calibration polynomial
            assert(restart==False)
            self.loadCalibration()
            self.justloaded = True
        except (IOError, AssertionError):
            self.speak('makine new calibration')
            self.justloaded = False
            self.create(**kw)
    
    def create(self, i=0, aperture=8.0, Teff=7000, N=50, remake=False, **kw):
        '''
        Populate the wavelength calibration for this aperture,
        using input from the user to try to link up lines between
        the measured arc spectra and the known line wavelengths.
        '''

        self.speak("populating wavelength calibration")
        
        # pick which exposure to use for the wavelength cal
        n = self.aperture.obs.fileprefixes['science'][i]
        
        # load that image
        self.aperture.loadExtracted(n)
        
        coef = find_wavecal(self.aperture.extracted, aperture=aperture, Teff=Teff, N=N, **kw)
        
        # convert nm to angstrom!
        self.pixelstowavelengths = np.poly1d(coef)

        self.saveCalibration()
        plt.show()
        
    def saveCalibration(self):
        '''Save the wavelength calibration.'''
        np.save(self.calibrationfilename, (self.pixelstowavelengths.coef))
        self.speak("saved wavelength calibration coefficients to {0}".format(self.calibrationfilename))

    def loadCalibration(self):
        '''Load the wavelength calibration.'''
        self.speak('trying to load calibration data')
        coef = np.load(self.calibrationfilename)
        self.pixelstowavelengths = np.poly1d(coef)
        self.polynomialdegree = len(coef) - 1
        self.speak("loaded wavelength calibration"
                "from {0}".format(self.calibrationfilename))