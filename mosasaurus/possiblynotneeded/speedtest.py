import TransmissionSpectrum
t = TransmissionSpectrum.TransmissionSpectrum('../wasp94_140805.obs')
b = t.bins[22]
b.fit( plot=False, slow=True, nburnin=500, ninference=1000)
'''def dofit(i):
    t.bins[i].fit(*t.psi(), plot=False, slow=True)

import multiprocessing
n = multiprocessing.cpu_count()
print "{n} cpus!".format(**locals())
pool = multiprocessing.Pool(n)
results = pool.map(dofit, range(4))'''


'''import cProfile, pstats
cProfile.run('b.fit(*t.psi(), plot=False, slow=True)', 'speedtest.profile')
p = pstats.Stats('speedtest.profile')
p.sort_stats('tottime').print_stats()
'''
