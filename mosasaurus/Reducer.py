from .imports import *
from .Observation import Observation
#from .CCD import CCD
from .Calibration import Calibration
from .Mask import Mask
from .Aperture import Aperture

class Reducer(Talker):
    '''Reducers are objects for reducing data, either interactively or not.

        This could be pictured as a mosasaurus,
        wearing a loupe and carrying lots of
        data structures around in its knapsack.'''

    def __init__(self, obs, label='default', visualize=True):
        '''initialize from an observation'''

        # decide whether or not this Reducer is chatty
        Talker.__init__(self)

        # should we visualize the steps of the reduction?
        self.visualize = visualize

        # connect to the observation
        self.obs = obs
        self.obs.reducer = self
        self.speak('the reducing mosasaurus is grabbing a fly spanker to analyze\n {0}'.format(self.obs))

        # setup all the other components of the reducer
        self.setup()
        self.speak('mosasaurus is ready to reduce')

        # set up the directory (inside the observation) to hold this extraction
        self.label=label
        self.extractionDirectory = os.path.join(self.obs.directory, "extraction_{}".format(self.label))
        mkdir(self.extractionDirectory)

    def setup(self):
        '''
        Setup all the basic components, or give them easier shortcuts.
        '''

        # give ourselves some shortcuts, in case we need them
        self.instrument = self.obs.instrument
        self.night = self.obs.night

        # create an interactive display to be associated with this
        self.display = loupe()

        # set up the Calibration
        self.calib = Calibration(self)

        # create a mask
        self.mask = Mask(self)

    def deleteandrestart(self):
        '''
        This deletes all intermediate data files (including log exposure choices)
        for this particular observation. Use this with great caution!
        '''

        self.speak("WARNING! You're about to erase all files in {0}".
                            format(self.obs.directory))

        if 'y' in self.input(" are you sure you're okay with this? [y,N]").lower():
            if 'y' in self.input("   for realsies? [y,N]").lower():
                shutil.rmtree(self.obs.directory)

    def reduce(self, remake=False):
        '''
        Process 2D multiobject spectral images into 1D spectra.
        '''

        filename = os.path.join(self.extractionDirectory, 'reductionstatus.txt')

        if os.path.exists(filename):
            self.speak('reduction in {} is already complete!'.format(self.extractionDirectory))
        else:
            self.speak('starting reductions for {}'.format(self.obs))

            # set up the calibrations and mask
            self.calib.setup()
            self.mask.setup()

            # loop over exposures and apertures
            self.mask.extractEverything(remake=remake)

            # write out that we're finished
            with open(filename, 'w') as f:
                f.write('Reduction was a success. Huzzah!')
