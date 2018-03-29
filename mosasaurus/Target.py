from .imports import *
from craftroom.star import Star

class Target(Talker):
    def __init__(self, starname='GJ1132', name='GJ1132b', **kw):
        '''
        Set up the target (at least the star with RA + Dec.)
        '''

        Talker.__init__(self)

        # what's the name of this target ("GJ1132b", planet, EB, etc...)
        self.name = name   

        # the starname might be different from the target (e.g. GJ1132)
        self.starname = starname

        # set up the star (to get its RA + Dec)
        self.setupStar(starname, **kw)

    def __repr__(self):
        '''How should this object be represented as a string?'''
        return '<Target {}>'.format(self.name)


    def setupStar(self, name, **kw):
        '''
        By default, use a name to initialize the star.
        Otherwise, pass ra + dec (etc...) keywords on to the Star object,
        to force it to initalize without querying SIMBAD.
        '''
        self.star = Star(name, **kw)

    def setupPlanet(self):
        '''Set up the initial conditions for planet fitting.'''
        pass
