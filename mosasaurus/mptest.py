import TransmissionSpectrum
import multiprocessing
obs = 'wasp94_140805.obs'
ncpu = multiprocessing.cpu_count()

def fastfit(i):
    t = TransmissionSpectrum.TransmissionSpectrum(obs)
    t.speak('starting fit for bin {0}'.format(i))
    t.bins[i].fit(plot=False, slow=False, interactive=False, remake=True)

def slowfit(i):
    t = TransmissionSpectrum.TransmissionSpectrum(obs)
    t.speak('starting fit for bin {0}'.format(i))
    t.bins[i].fit(plot=False, slow=True, interactive=False, nburnin=500, ninference=500)

pool = multiprocessing.Pool(ncpu)
t = TransmissionSpectrum.TransmissionSpectrum(obs)
for i in range(len(t.bins)):
    fastfit(i)
#pool.map_async(fastfit, range(len(t.bins)))
#pool.map_async(slowfit, range(len(t.bins)))
