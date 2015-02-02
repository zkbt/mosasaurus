from imports import *

class Mask():
  '''Mask object keeps track of the collection of apertures that make up the mask.'''
  def __init__(self, calib):
    '''Initialize the mask object.'''
    self.calib = calib
    self.ccd = self.calib.ccd
    self.obs = self.calib.obs
    self.display = self.obs.display
    self.setup()

  def extractCenters(self, ds9RegionString):
    '''Extract x and y pixel positions from a string of regions returned by pyds9.'''
    regions = ds9RegionString.split('circle')[1:]
    x, y = [],[]
    for region in regions:
      x.append(float(region[1:].split(',')[0]))
      y.append(float(region[1:].split(',')[1]))
    return np.array(x), np.array(y)

  def pickStars(self):
    '''Open ds9 and pick the stars we want to reduce.'''
    extractionCentersFilename = self.obs.workingDirectory + 'extractionCenters.txt'
    if os.path.exists(extractionCentersFilename):
      print "Looks like extraction centers were already defined; loading them from " + extractionCentersFilename
      x, y = np.loadtxt(extractionCentersFilename)
    else:
      print "You should pick the stars you're interested in."
      self.obs.display.one(self.calib.images['Undispersed'], clobber=True)
      D = self.obs.display.ds9
      D.set("scale log")
      D.set("regions centroid auto yes")
      var = raw_input("Plop down regions on the stars you'd like to extract. Then hit enter:\n  (ds9 should auto-centroid)\n")
      regions = D.get("regions")
      x,y = self.extractCenters(regions)
      np.savetxt(extractionCentersFilename, (x,y))
    self.xextract, self.yextract = x, y
    return x, y


  def setup(self, visualize=True, remake=False):
    '''Setup extraction masks, starting with the x and y positions selected by user.'''

    # open tool for use to pick stars
    self.pickStars()

    try:
      # ultimately, include ability to load extraction parameters that were already defined
      assert(False)
    except:
      #load the apertures, if they exist already
      if os.path.exists(self.obs.workingDirectory + 'apertures.npy') and remake==False:
        apertures = np.load(self.obs.workingDirectory + 'apertures.npy')
        return apertures

      self.apertures = []
      # open a dispersed image (take the first one)
      imageDispersed = self.calib.science()


      for i in range(len(self.xextract)):
        x, y = self.xextract[i], self.yextract[i]
        a = Aperture(x,y,self)
        self.apertures.append(a)

      if visualize:
        r = zachopy.regions.Regions('aperture')
        for a in self.apertures:
          # display the target positions on the image
          r.addCircle(a.x, a.y, text="({0:.0f}, {1:.0f})".format(x, y), color='green', font="helvetica 10 bold")
          r.addBox(a.xbox, a.ybox, a.wbox, a.hbox, color='green')
        filename = self.obs.extractionDirectory + 'apertures.reg'
        r.write(filename)
        self.obs.display.one(imageDispersed, clobber=True)
        self.obs.display.ds9.set("regions load {0}".format(filename))


  def createStamps(self, n):
    '''Write postage stamps around the spectra.'''
    print "   creating (unfiltered) postage stamps for all science images"
    #for n in self.obs.nScience:
    print "     cutting stamps out of ccd{0:04}".format(n)
    this = self.ccd.readData(n)
    for a in self.apertures:
      stamp = a.stamp(this)
      writeFitsData(stamp, stampFilename(n,a))
      print "       saved stamp to", stampFilename(n, a)
