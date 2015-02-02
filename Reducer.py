# required things to add:
# estimate gain for amplifiers, convert everything to electrons early on

from imports import *
from Observation import Observation
from Night import Night
from Tools import *
from Display import Display
from CCD import CCD
from Calibration import Calibration
from Mask import Mask
from Aperture import Aperture

class Reducer(Talker):
    '''Object for reducing data, either interactively or not. Picture as a mosasaurus,
        wearing a loup and carrying lots of data structures around in its knapsack.'''

    def __init__(self, filename='wasp94_140805.obs', **kwargs):
        # decide whether or not this Reducer is chatty
        Talker.__init__(self, **kwargs)

        self.speak('a new data-reduction mosasaurus has been created')

        # load in the observation file for this object
        self.filename = filename
        self.obs = Observation(self.filename)

        self.speak('mosasaurus is ready to reduce')

    def reduce(self):
        self.speak('reducing observations of {0} on {1}'.format(self.obs.name, self.obs.night))

        # load the headers for this observation
        self.obs.loadHeaders()

        # create a night object, associated with this observation
        night = Night(self.obs)
        # create an observing log from all the image headers from this night
        night.obsLog()

        # set up the calibration
        calib = Calibration(self.obs)

        '''# loop through the CCD's needed for this observation, and make sure they are stitched
        #for n in obs.nNeeded:
        #		ccd = CCD(obs,n=n,calib=calib)
        #	ccd.createStitched(visualize=True)


        mask = Mask(calib)
        for a in mask.apertures:
            a.displayStamps(a.images)
            a.extractAll(remake=True)'''
