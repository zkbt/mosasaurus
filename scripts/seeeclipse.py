import matplotlib.pyplot as plt
import zachopy.oned, numpy as np
trans = np.load('onetrans.npy')[()]
plt.ion()
bins = {}
albedo=1
def scale(y):
    return (y-1)*1e6
f = plt.figure('update')
axupdate = plt.subplot()
e = plt.figure('bin')
axbin = plt.subplot()
count=0
for b in trans.bins:
    w, t, residuals, ok = b.wavelength,\
        b.tlc.timefrommidtransit(),\
        b.tlc.residuals(),\
        b.tlc.bad == False


    duration = b.tm.planet.duration
    center = 0.0 + duration*2./3
    model = np.ones_like(t)
    model[np.abs(t - center) < duration/2.0] = 1.0 - 220.0*1e-6*albedo
    data = model + residuals
    std = np.std(residuals)

    '''try:
        stack, weights
    except:
        stack, weights = np.zeros_like(t), np.zeros_like(t)
    stack[ok] += data[ok]/std**2
    weights[ok] += np.ones_like(data)[ok]/std**2
    axupdate.cla()
    axupdate.plot(t[ok], stack[ok]/weights[ok])'''
    plt.cla()

    axbin.plot(t[ok], scale(data[ok]), color='gray', alpha=0.5)
    axbin.plot(t, scale(model), color='red', alpha=0.5)
    bx, by, be = zachopy.oned.binto(t[ok], data[ok], binwidth=0.01)
    axbin.plot(bx, scale(by), marker='o', alpha=0.5, color='black', linewidth=0)
    ax = plt.gca()
    ax.set_ylim(-5000,5000)
    bins[w] = (t,residuals)

    try:
        array
    except:
        array = np.zeros((len(by), len(trans.bins)))
        n = array.shape[0]-1
    array[:n,count] = by[:n]
    print array

    a = raw_input(str(w))
    count += 1
