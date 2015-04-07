import mosasaurus.transmission as transmission
import matplotlib.pyplot as plt, numpy as np

files = ['wasp94_140801.obs', 'wasp94_140805.obs', 'wasp94_140809.obs']
spectra = []
plt.figure('three spectra')
toshow = ['k']
plt.figure('comparison')
gs = plt.matplotlib.gridspec.GridSpec(1,1)
count =0
sharex=None
ax = {}
medians = {}
for file in files:
    t = transmission.load(file, 100, 'fixedGeometry')
    t.load()
    spectra.append(t)
    w = np.array([b.wavelength for b in t.bins])/t.unit
    for k in toshow:
        try:
            ax[k]
            medians[k]
        except:
            ax[k] = plt.subplot(gs[count], sharex=sharex)
            count += 1
            ax[k].set_ylabel(k)
            sharex = ax[k]
            medians[k] = []


        p = np.array([b.tm.planet.__dict__[k].value for b in t.bins])
        u = np.array([b.tm.planet.__dict__[k].uncertainty for b in t.bins])
        medians[k].extend(p)
        print k
        print np.median(medians[k])
        ax[k].errorbar(w, p, u, marker='o', markersize=10, linewidth=3, elinewidth=3, capsize=5, capthick=3, alpha=0.5)

    plt.draw()
    a = raw_input('!!!')
