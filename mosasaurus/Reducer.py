# required things to add:
# estimate gain for amplifiers, convert everything to electrons early on

from imports import *
from Observation import Observation
from Night import Night
from Tools import *
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

    def reduce(self, remake=False):
        self.speak('reducing observations of {0} on {1}'.format(self.obs.name, self.obs.night))

        # load the headers for this observation
        self.obs.loadHeaders()

        # create a night object, associated with this observation
        night = Night(self.obs)

        # create an observing log from all the image headers from this night
        night.obsLog()

        # set up the Calibration
        self.speak('setting up Calibration data')
        self.calib = Calibration(self.obs)

        # loop through the CCD's needed for this observation, and make sure they are stitched
        self.speak('stitching all science images and rejecting cosmics along the way')
        #self.calib.rejectCosmicRays()

        # create a mask
        self.mask = Mask(self.calib)

        # loop over exposures and apertures
        self.mask.extractEverything(remake=remake)
