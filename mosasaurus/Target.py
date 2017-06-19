from imports import *
from zachopy.star import Star

class Target(Talker):
    def __init__(self, starname='GJ1132', name='GJ1132b', **kw):
        '''
        Set up the target (at least the star with RA + Dec.)
        '''

        Talker.__init__(self)

        # what's the name of this target (planet, EB, etc...)
        self.name = name

        # set up the star (to get its RA + Dec)
        self.setupStar(starname, **kw)

    def setupStar(self, name, **kw):
        '''
        By default, use a name to initialize the star.
        Otherwise, pass ra + dec (etc...) keywords on to the Star object.
        '''
        self.star = Star(name)

    def setupPlanet(self):
        '''Set up the initial conditions for planet fitting.'''
        pass
