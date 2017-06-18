import TransmissionSpectrum
import multiprocessing

kw = dict(obs = '../wasp94_140801.obs', label='floatingGeometry')
ncpu = multiprocessing.cpu_count()


# create a base transmission spectrum object
t = TransmissionSpectrum.TransmissionSpectrum(**kw)
t.loadLCs()
n = len(t.bins)
inputs = [(i,kw) for i in range(n)]
print inputs[0]
# do a fast fit, to figure out where the outliers are
map(TransmissionSpectrum.fastfit, inputs)


# mask the outliers
t.load(method='lm')
t.createMask(afterfastfit=True)

kw['maskname'] = 'trimOutliers'
inputs = [(i,kw) for i in range(n)]
pool = multiprocessing.Pool(ncpu)
pool.map_async(TransmissionSpectrum.slowfit, inputs)
#pool.close()
#pool.join()
