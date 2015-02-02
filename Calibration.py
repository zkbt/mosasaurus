from imports import *
from CCD import CCD
plt.ion()
class Calibration(Talker):
  '''Calibration object stores information related to the calibration for this observation.'''
  def __init__(self, obs, visualize=True, **kwargs):
    '''Initialize calibration object.'''

    # decide whether or not this Reducer is chatty
    Talker.__init__(self, **kwargs)

    self.obs = obs
    self.display = zachopy.display.ds9('calibration')
    self.ccd = CCD(self.obs, calib=self)
    self.visualize = visualize


    self.createBadPixelMask()
    self.createMasterImages()

  def estimateGain(self):
    try:
      self.gains = np.loadtxt(self.obs.workingDirectory + 'gains.txt')
      self.speak("loaded gains from {0}".format(self.obs.workingDirectory + 'gains.txt'))
      self.speak("   they are:")
      for g in self.gains:
        self.speak("       {0}".format(g))
    except:
      self.speak("estimating gain from noise in multiple flat-field exposures.")
      c = CCD(self.obs, calib=self)

      fi = plt.figure('gain estimation', figsize=(7,10))
      gs =plt.matplotlib.gridspec.GridSpec(1,2, hspace=0, wspace=0)
      ax = []
      ax.append(plt.subplot(gs[0]))
      ax.append(plt.subplot(gs[1], sharey=ax[0], sharex=ax[0]))
      plt.setp(ax[1].get_yticklabels(), visible=False)
      s1, s2 = [], []
      for n in self.obs.nWideFlat:
        c.set(n, 'FlatInADU')
        c.readData()
        c1, c2 = c.amplifiers()
        s1.append(c1)
        s2.append(c2)
      s1 = np.array(s1)
      s2 = np.array(s2)
      gains = [1.0, 1.0]
      for i in [0,1]:
        s = [s1, s2][i]
        name = ['c1', 'c2'][i]
        self.speak("")
        self.speak(name)
        self.speak("")
        medianinadu = np.median(s, 0).flatten()
        #noiseinadu = 1.48*np.median(np.abs(s - medianinadu.reshape(1,s.shape[1],s.shape[2])), 0).flatten()
        noiseinadu = np.std(s[1:,:,:]-s[0:-1,:,:], 0).flatten()/np.sqrt(2)
        ok = (medianinadu > 25000)*(medianinadu < 45000)

        def noisemodel(n, gain):
          return np.sqrt(n)/np.sqrt(gain)

        def deviates(parameters, n=None, noise=None, fjac=None):
          gain = parameters[0]
          status = 0
          dev = np.log(noise/noisemodel(n,gain))
          #dev = (noise/noisemodel(n, gain) - 1)
          noiseondev = zachopy.utils.mad(dev)
          normalized = dev/noiseondev
          return [status, normalized]



        converged = False
        oldok = ok + 0
        while(converged == False):
          p0 = [1.0]
          parinfo = [{'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]} for n in range(len(p0))]
          parinfo[0]['limits']=[0.1, 10.0]
          parinfo[0]['limited']=[1,1]
          n = medianinadu
          noise = noiseinadu
          fit = mpfit(deviates, p0, parinfo=parinfo, functkw={'n':n[ok], 'noise':noise[ok]})

          gain = fit.params[0]



          plt.draw()
          ok *= (np.abs(deviates(fit.params, n=n, noise=noise)[-1]) < 5)
          converged = (ok == oldok).all()
          oldok = ok + 0
          self.speak("{0} bad points! converged? {1}".format(np.sum(ok == False), converged))
          self.speak( "   {0}".format(gain))
        gains[i] = gain

        ax[i].cla()
        ax[i].plot(n[ok].flatten(),noise[ok].flatten(), marker='o', linewidth=0, alpha=0.1, color='black', markersize=3)
        ax[i].plot(n[ok == False].flatten(),noise[ok == False].flatten(), marker='o', linewidth=0, alpha=0.02, color='red', markersize=3, markeredgewidth=0)
        x = np.linspace(0, n[ok].max(), 100)
        ax[i].plot(x, noisemodel(x,gain), color='green', linewidth=5, alpha=0.5)
        ax[i].plot(x, noisemodel(x,1.5), color='blue', linewidth=5, alpha=0.2)
        ax[i].plot(x, noisemodel(x,1.8), color='blue', linewidth=5, alpha=0.2)
        ax[i].set_title('c{0} gain estimated to be {1:.2f}'.format(i, gain))



        #ax[i].plot(n[ok].flatten(), deviates(fit.params, n=n[ok], noise=noise[ok])[-1].flatten(), marker='o', linewidth=0, alpha=0.005, color='black', markersize=3)
        #ax[i].plot(n[ok==False].flatten(), deviates(fit.params, n=n[ok == False], noise=noise[ok == False])[-1].flatten(), marker='o', linewidth=0, alpha=0.1, color='red', markersize=3)
        self.speak('Does this gain estimate seem reasonable?')
        a = raw_input(name)
      np.savetxt(self.obs.workingDirectory + 'gains.txt', gains)
      self.gains = gains

  def createMedianImage(self, n, visualize=True, type=None):
    '''Take median of a series of images (requires enough memory to hold them all.)'''
    array = self.ccd.loadImages(n, type=type)
    if len(array.shape) <=2:
      return array
    median = np.median(array,0)
    #if visualize:
      #self.display.many(array, clobber=True)
      #self.display.one(median, clobber=False)
    return median

  def createMeanImage(self, n, cosmic=True, visualize=False, type=None):
    '''Take the mean of a series of images, less memory intensive than median.'''
    stride = np.maximum(len(n)/100, 1)
    for i in np.arange(0, len(n), stride):
      data = self.ccd.readData(n[i], type=type)
      if cosmic and i == 0:
        outlier = np.zeros_like(data)
      if i == 0:
        count = 1
        summedImage = data
        if cosmic:
          summedSquaredImage = data**2
          stddev = np.sqrt(np.maximum(summedSquaredImage/count - (summedImage/count)**2,1))
      else:
        count += 1
        summedImage += data
        if cosmic:
          summedSquaredImage += data**2
          stddev = np.sqrt(np.maximum(summedSquaredImage/count - (summedImage/count)**2,1))
          bad = (data - last)/stddev > self.obs.cosmicThreshold
          outlier[bad] += (data[bad] - summedImage[bad]/count)
      last = data
      self.speak('        ' + self.ccd.name)
      #self.display.one(summedImage, clobber=(i == 0))
      if cosmic:
        #self.display.one((data - summedImage/count)/stddev)
        if visualize:
          if i == 0:
            self.display.one(outlier, clobber=True)
          else:
            self.display.ds9update(outlier)
    if visualize:
      self.display.one(summedImage)
      self.display.one(summedImage - outlier)

    if cosmic:
      cosmicFilename = self.obs.workingDirectory + 'cosmics{0}to{1}.fits'.format(np.min(n), np.max(n))
      writeFitsData(outlier, cosmicFilename)
    return (summedImage-outlier)/count, stddev



  def createBadPixelMask(self, visualize=True):
    '''Try to estimate bad pixels from a flat image. KLUDGE'''
    self.speak("populating bad pixel mask")
    badPixelFilename = self.obs.workingDirectory + 'badpixels.fits'
    try:
      self.images['BadPixels'] = readFitsData(badPixelFilename)
      self.speak( "loaded bad pixel mask from ", badPixelFilename)
    except:
      self.speak( "creating bad pixel mask from the master flat frames")
      c = CCD(self.obs, calib=self)

      cube = []
      for n in self.obs.nWideFlat:
        c.set(n, 'WideFlat')
        cube.append(c.readData())

      cube = np.array(cube)

      median = np.median(cube,0)
      noise = np.median(np.abs(cube - median.reshape(1,cube.shape[1], cube.shape[2])), 0)
      plt.figure('bad pixel mask')
      ax = plt.subplot()
      ax.plot(median.flatten(), noise.flatten(), color='black', alpha=0.5, marker='o', markersize=4, markeredgewidth=0, linewidth=0)
      ax.set_yscale('log')
      bad = (noise < 0.05*np.sqrt(median)) | (noise == 0) | (median == 0) | (median < 0) | (self.bias() > 10000) | (self.dark() > 100)
      ax.plot(median[bad].flatten(), noise[bad].flatten(), color='red', alpha=0.5, marker='o', markersize=10, markeredgecolor='red', linewidth=0)

      self.images['BadPixels'] = bad.astype(np.int)
      if visualize:
        self.display.one(self.images['BadPixels'])
        self.input("Do you approve of the bad pixel mask?")
      writeFitsData(self.images['BadPixels'], badPixelFilename)

  def createMasterImage(self, type=None, remake=False):
    try:
      self.images.keys()
    except:
      self.images = {}

    # if no name included, do nothinh
    if type is None:
      return

    self.ccd.set(n=None, type=type)
    # make sure this is a master image that is possible to make
    assert(type in self.obs.cal_dictionary.keys())

    self.speak( "   populating the master {0} image".format(type))
    masterFilename = self.obs.workingDirectory + "master_{0}.fits".format(type)
    try:
      if type == 'Science':
        self.images['Science'] = readFitsData(self.obs.workingDirectory + "master_{0}.fits".format('Science'))
        self.images['ScienceStdDev'] = readFitsData(self.obs.workingDirectory + "master_{0}.fits".format('ScienceStdDev'))
      else:
        self.images[type] = readFitsData(masterFilename)
      self.speak( "     loaded from {0}".format( masterFilename))
    except:
      self.speak("     creating from images " + zachopy.utils.truncate( str( self.obs.cal_dictionary[type]), n=30))
      if type == 'Science':
        self.images['Science'], self.images['ScienceStdDev'] = self.createMeanImage(self.obs.cal_dictionary[type], type=type)
        writeFitsData(self.images['Science'], self.obs.workingDirectory + "master_{0}.fits".format('Science'))
        writeFitsData(self.images['ScienceStdDev'], self.obs.workingDirectory + "master_{0}.fits".format('ScienceStdDev'))
      else:
        self.images[type] = self.createMedianImage(self.obs.cal_dictionary[type], type=type)
        writeFitsData(self.images[type], masterFilename)

    self.display.one(self.images[type])
    self.input("Do you like master image {0}?".format(type))

  def bias(self):
    try:
      return self.images['Bias']
    except:
      self.createMasterImage('Bias')
      return self.images['Bias']

  def dark(self):
    try:
      return self.images['Dark']
    except:
      self.createMasterImage('Dark')
      return self.images['Dark']

  def science(self):
    try:
      return self.images['Science']
    except:
      self.createMasterImage('Science')
      return self.images['Science']

  def wideflat(self):
    try:
      return self.images['WideFlat']
    except:
      self.createMasterImage('WideFlat')
      return self.images['WideFlat']


  def createMasterImages(self, remake=False):
    '''Combine individual exposures into master frames for the various calibrations.'''

    for k in self.obs.cal_dictionary.keys():
      self.createMasterImage(k, remake=remake)
