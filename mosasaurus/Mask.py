from .imports import *
from .Aperture import Aperture

class Mask(Talker):
  '''Mask objects keep track of the collection of Apertures
        that should be looked at for this observation.

            This can be thought of as the mosasaurus' checklist
            of slits to investigate, and bookkeeping tools.'''

  def __init__(self, reducer, **kwargs):
    '''Initialize the mask object.'''

    # decide whether or not this Mask is chatty
    Talker.__init__(self, **kwargs)

    # connect it to other parts
    self.reducer = reducer
    self.calib = self.reducer.calib
    self.ccd = self.calib.ccd
    self.obs = self.calib.obs

    # set up a display
    self.display = self.reducer.display
    self.speak('created a mask, to store extraction regions')
    #self.setup()

  def pickStars(self):
    '''Open ds9 and pick the stars we want to reduce.'''
    extractionCentersFilename = os.path.join(self.reducer.extractionDirectory, 'extractionCenters.txt')
    if os.path.exists(extractionCentersFilename):
      self.speak("Extraction centers were already defined; loading them from {0}".format(extractionCentersFilename))
      x, y = np.transpose(np.loadtxt(extractionCentersFilename))
    else:
      self.speak("Please should pick the stars you're interested in, then [q]uit.")

      # create empty list of extraction centers
      self.xcenters, self.ycenters = [], []

      # set up a loupe display to show an Reference image
      self.loupe = self.reducer.display

      # add in an option to pull out an extraction center
      self.loupe.options['a'] = dict(description='[a]dd an extraction center',
                                        function=self.addExtractionCenter,
                                        requiresposition=True)

      # add in an option to pull out an extraction center
      self.loupe.options['r'] = dict(description='[r]estart, clearing all centers',
                                        function=self.resetExtractionCenters,
                                        requiresposition=True)

      # set up the plots
      self.loupe.setup(self.calib.images['reference'])

      # start up the loupe event handler
      self.loupe.run()

      # pull out the extraction centers and save them
      y,x = self.xcenters, self.ycenters
      assert(len(x) > 0)
      np.savetxt(extractionCentersFilename, np.transpose([x,y]))

    self.xextract, self.yextract = x, y
    self.speak('extraction centers are')
    self.speak('   x={}'.format(x))
    self.speak('   y={}'.format(y))


    return x, y


  def resetExtractionCenters(self, pressed):
      self.xcenters, self.ycenters = [], []
      try:
        self.centersplot.set_data(self.xcenters, self.ycenters)
      except:
        self.centersplot = self.loupe.ax['2d'].plot(self.xcenters, self.ycenters,
                                                    linewidth=0, marker='x',
                                                    markersize=10, color='black')

  def addExtractionCenter(self, pressed):
      '''from a keyboard event, add an extraction center'''
      x, y = pressed.xdata, pressed.ydata
      self.speak('adding an extraction center at {:.1f}, {:.1f}'.format(x,y))
      self.xcenters.append(x.astype(np.int))
      self.ycenters.append(y.astype(np.int))
      self.speak(' {} {}'.format(self.xcenters, self.ycenters))

      try:
        self.centersplot.set_data(self.xcenters, self.ycenters)
      except:
        self.centersplot = self.loupe.ax['2d'].plot(self.xcenters, self.ycenters,
                                                    linewidth=0, marker='x',
                                                    markersize=10, color='black')



  def setup(self, visualize=True, remake=False):
    '''Setup extraction masks, starting with the x and y positions selected by user.'''

    # open tool for use to pick stars
    self.pickStars()

    try:
      # ultimately, include ability to load extraction parameters that were already defined
      assert(False)
    except:
      #load the apertures, if they exist already
      if os.path.exists(self.reducer.extractionDirectory + 'apertures.npy') and remake==False:
        apertures = np.load(self.reducer.extractionDirectory + 'apertures.npy')
        return apertures

      self.apertures = []
      # open a dispersed image (take the first one)
      imageDispersed = self.calib.science()


      for i in range(len(self.xextract)):
        x, y = self.xextract[i], self.yextract[i]
        a = Aperture(x,y,self)
        self.apertures.append(a)

      filename = os.path.join(self.reducer.extractionDirectory, 'genericfinderchart.pdf')
      if visualize and (os.path.exists(filename) == False):
        # create a finder chart of apertures
        self.loupe = self.reducer.display
        self.loupe.one(self.calib.images['reference'])

        for a in self.apertures:
          # display the target positions on the image
          self.loupe.ax['2d'].plot(a.y, a.x, marker='x', color='black', alpha=0.5, markersize=8, linewidth=0)
          self.loupe.ax['2d'].text(a.y, a.x, '   {:03.0f},{:03.0f}'.format(a.x, a.y), color='black', alpha=0.5, ha='left', va='center', fontsize=6)


          #r.addBox(a.xbox, a.ybox, a.wbox, a.hbox, color='green')
        plt.savefig(filename)
        self.speak('saved finder chart to {}'.format(filename))


  def createStamps(self, n):
    '''Write postage stamps around the spectra.'''
    self.speak("   creating (unfiltered) postage stamps for all science images")
    #for n in self.obs.nScience:
    self.speak("     cutting stamps out of ccd{0:04}".format(n))
    this = self.ccd.readData(n)
    for a in self.apertures:
      stamp = a.stamp(this)
      writeFitsData(stamp, stampFilename(n,a))
      self.speak("saved stamp to {0}".format(stampFilename(n, a)))


  def load(self, n):
      '''Make sure the mask's CCD data is set to the correct frame.'''
      self.ccd.readData(n, imageType='science')
      assert(self.ccd.exposureprefix == n)
      self.speak('set CCD data to {0}'.format(self.ccd.name))


  def extractStars(self, n, remake=False):
      '''For one frame, extract all the spectra in this mask.'''

      # provide update
      self.speak('extracting {0} spectra from image #{1}'.format(len(self.apertures), n))

      # load the (entire!) image
      #self.load(n)

      # loop over apertures
      for a in self.apertures:
          # extract the spectrum in this aperture
          a.extract(n, remake=remake)

      # hzdl - testing:
      #self.apertures[4].extract(n, remake=remake)

          #self.input('just finished extracting all stars ({0})'.format(n))

      #import sys
      #sys.exit("Breaking here. Check it out.")

  def extractEverything(self, remake=False):

    '''Loop through exposures, extracting all the spectra in this mask.'''
    self.speak('looping through all frames and extracting all spectra in them')

    # extract spectra
    self.speak('making sure all spectra have been extracted')
    for exposureprefix in self.obs.fileprefixes['science']:
        self.extractStars(exposureprefix, remake=remake)

    # make a movie of the extractions
    if False:
        self.speak('making movies of the extraction apertures')
        for a in self.apertures:
            a.movieExtraction()


    # add wavelength calibration
    self.speak('adding wavelength calibrations to all spectra (if not already done)')
    self.addWavelengthCalibration()

  def addWavelengthCalibration(self, remake=False, shift=False):
      '''don't extract, just addWavelengthCalibration the wavelengths and resample'''

      self.speak('creating wavelength calibrators')
      for aperture in self.apertures:
          aperture.createWavelengthCal(remake=remake)

      self.speak('adding wavelength calibrations to all stars (if needed)')
      for exposureprefix in self.obs.fileprefixes['science']:
          for a in self.apertures:
              # a.visualize = False
              a.addWavelengthCalibration(exposureprefix, shift=shift)
