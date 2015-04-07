import transmission
import matplotlib.pyplot as plt, numpy as np
from astropy.io import ascii
files = ['wasp94_140801.obs', 'wasp94_140805.obs', 'wasp94_140809.obs']
toshow = ['k']

# an empty dictionary of observations
spectra = {}

# an empty dictionary of bins (which will each contain a list of spectra)
bins = {}

# create a figure
plt.figure('three spectra')
gs = plt.matplotlib.gridspec.GridSpec(1,1)
sharex=None
ax = {}
count =0


fort = ascii.read("lambda_1500K_g10_wTiOVO.dat")
fort_lam = fort['wavelength']*10000
fort_radius= fort['radiusinkm']
fort_radius = fort_radius[fort_lam.argsort()]
fort_lam.sort()
r_jupiter = 7149200000 # cm
r_sun = 69550000000 # cm
fort_rp_over_rs = fort_radius*100000.0*1.53/1.37/(1.44*r_sun)
fort_depth = 100*fort_rp_over_rs**2
fort_depth = fort_depth* 1.2/np.mean(fort_depth[(fort_lam >7000)*(fort_lam<10000)])

for file in files:

    # load the transmission spectrum
    t = transmission.load(file, 100, 'fixedGeometry')
    t.load()

    spectra[file] = t

    w = np.array([b.wavelength for b in t.bins])/t.unit
    for k in toshow:
        try:
            ax[k]
        except:
            ax[k] = plt.subplot(gs[count], sharex=sharex)
            count += 1
            ax[k].set_ylabel(k)
            sharex = ax[k]


        p = np.array([b.tm.planet.__dict__[k].value for b in t.bins])
        u = np.array([b.tm.planet.__dict__[k].uncertainty for b in t.bins])
        print file
        print k
        print 2*p*u*1e6
        print 'LD'
        print ([b.tm.star.__dict__['u1'].value for b in t.bins])
        ax[k].errorbar(w, p, u, marker='o', markersize=10, linewidth=3, elinewidth=3, capsize=5, capthick=3, alpha=0.25)
    for i in range(len(w)):
        this = {k:p[i], k+'_uncertainty':u[i]}
        try:
            bins[w[i]].append(this)
        except:
            bins[w[i]] = [this]



    plt.draw()
    a = raw_input('!!!')

binned, unc = {}, {}
for k in bins.keys():
    values = np.array([b['k'] for b in bins[k]])
    uncertainties = np.array([b['k_uncertainty'] for b in bins[k]])
    binned[k] = np.sum(values/uncertainties**2)/np.sum(1.0/uncertainties**2)
    unc[k] = np.sqrt(1/np.sum(1.0/uncertainties**2))

k = binned.keys()
ax = plt.gca()
ax.errorbar(np.array(k), np.array([binned[i] for i in k]), np.array([unc[i] for i in k]),  marker='o', markersize=10, linewidth=0, elinewidth=3, capsize=5, capthick=3, alpha=0.5, color='black')
ax = plt.gca()
ax.plot(fort_lam/10.0, np.sqrt(fort_depth/100.0), zorder=-100, alpha=0.25, color='gray')
