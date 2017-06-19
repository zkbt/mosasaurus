# create a target
t = Target(starname='GJ1132', name='GJ1132b')

# create an instrument
i = LDSS3C(grism='vph-red')

# create object to identify data and organize this night
n = Night('ut160227_28', instrument=i)

# create an observation of this target, with this instrument, on this night
o = Observation(t, i, n)
