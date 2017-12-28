from .imports import *
from Observation import Observation
from Night import Night
from Tools import *
from CCD import CCD
from Calibration import Calibration
from Mask import Mask
from Aperture import Aperture

class Reducer(Talker):
    '''Reducers are objects for reducing data, either interactively or not.

        This could be pictured as a mosasaurus,
        wearing a loup and carrying lots of
        data structures around in its knapsack.'''

    def __init__(self, filename='wasp94_140805.obs', **kwargs):
        '''initialize from a ".obs" file'''

        # decide whether or not this Reducer is chatty
        Talker.__init__(self, **kwargs)

        self.speak('the reducing mosasaurus is grabbing a fly swatter to analyze {0}'.format(filename))
        # store the filename of the reduction parameter file
        self.filename = filename

        # setup all the components of the reducer
        self.setup()

        self.speak('mosasaurus is ready to reduce')

    def setup(self):

        # load in the observation file for this object
        self.obs = Observation(self.filename)

        self.display = loupe() #''mosasaurus',
                               # xsize=self.obs.xsize*self.obs.displayscale,
                               # ysize=self.obs.ysize*self.obs.displayscale)

        # create a night object, associated with this observation
        self.night = Night(self.obs)

        # set up the Calibration
        self.calib = Calibration(self)

        # create a mask
        self.mask = Mask(self)

    def deleteandrestart(self):
        self.speak("WARNING! You're about to erase all files in {0}".
                            format(self.obs.workingDirectory))

        if 'y' in self.input(" are you sure you're okay with this? [y,N]").lower():
            if 'y' in self.input("   for realsies? [y,N]").lower():
                shutil.rmtree(self.obs.workingDirectory)

        self.obs = Observation(self.filename)

    def reduce(self, remake=False):
        '''process 2D multiobject spectral images into 1D spectra'''

        self.speak('reducing observations of {0} on {1}'.format(
                        self.obs.name, self.obs.night))

        # set up the calibrations and mask
        self.calib.setup()
        self.mask.setup()

        # load the headers for this observation
        self.obs.loadHeaders()

        # create an observing log from all the image headers from this night
        self.night.obsLog()

        # loop over exposures and apertures
        self.mask.extractEverything(remake=remake)
