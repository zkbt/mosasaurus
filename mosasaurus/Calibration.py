from imports import *
from CCD import CCD

class Calibration(Talker):
  '''Calibrations are objects that store calibration data,
        including both afternoon exposures (biases, darks, flats)
        and some on-sky exposures (direct images, master spectral images).

            This can be thought of as bookmarked pages of
            the reducing mosasaurus' reference books,
            keeping track of calibrations that might be useful'''
  def __init__(self, reducer, visualize=True, **kwargs):
    '''Initialize calibration object.'''

    # decide whether or not this Reducer is chatty
    Talker.__init__(self, **kwargs)


    self.speak('setting up calibrator')
    self.reducer = reducer
    self.obs = self.reducer.obs
    self.display = self.reducer.display
    self.ccd = CCD(self.obs, calib=self)
    self.visualize = visualize
    self.images = {}

  def setup(self):
    self.createBadPixelMask()
    self.createMasterImages()
    self.speak('calibration data are processed and ready for use')

  def estimateGain(self):
    self.gains = self.obs.gains
    if self.gains is None:
        try:
          self.gains = np.loadtxt(self.obs.workingDirectory + 'gains.txt')
          self.speak("loaded gains from {0}".format(self.obs.workingDirectory + 'gains.txt'))
          self.speak("   they are:")
          for g in self.gains:
            self.speak("       {0}".format(g))
        except:
          self.speak("estimating gain from noise in multiple flat-field exposures.")
          c = self.ccd#CCD(self.obs, calib=self)

          fi = plt.figure('gain estimation', figsize=(10,4))
          gs =plt.matplotlib.gridspec.GridSpec(1,2, hspace=0, wspace=0, top=0.85)
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
          readnoises = [10.0, 10.0]
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
            ok = (medianinadu > 20000)*(medianinadu < 30000)

            def noisemodel(n, gain, readnoise, lamp):
              return np.sqrt(n/gain + readnoise**2 +lamp**2*n**2)

            def deviates(parameters, n=None, noise=None, fjac=None):
              gain = parameters[0]
              readnoise = parameters[1]
              lamp = parameters[2]
              status = 0
              dev = np.log(noise/noisemodel(n,gain,readnoise, lamp))
              #dev = (noise/noisemodel(n, gain,readnoise) - 1)
              noiseondev = zachopy.utils.mad(dev)
              normalized = dev/noiseondev
              return [status, normalized]



            converged = False
            oldok = ok + 0
            while(converged == False):
              p0 = [1.5, 10.0, 0.000]
              parinfo = [{'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]} for n in range(len(p0))]
              parinfo[0]['limits']=[0.1, 10.0]
              parinfo[0]['limited']=[1,1]

              parinfo[1]['limits']=[5.0, 20.0]
              parinfo[1]['limited']=[1,1]
              parinfo[1]['fixed']=1

              parinfo[2]['limits']=[0.0, 0.2]
              parinfo[2]['limited']=[1,1]
              parinfo[2]['fixed']=1

              n = medianinadu
              noise = noiseinadu
              fit = mpfit(deviates, p0, parinfo=parinfo, functkw={'n':n[ok], 'noise':noise[ok]})

              gain = fit.params[0]
              readnoise = fit.params[1]
              lamp = fit.params[2]


              plt.draw()
              ok *= (np.abs(deviates(fit.params, n=n, noise=noise)[-1]) < 5)
              converged = (ok == oldok).all()
              oldok = ok + 0
              self.speak("{0} bad points! converged? {1}".format(np.sum(ok == False), converged))
              self.speak( "   {0}".format(gain))
            gains[i] = gain
            readnoises[i] = readnoise
            ax[i].cla()
            ax[i].plot(n[ok].flatten(),noise[ok].flatten(), marker='o', linewidth=0, alpha=0.1, markeredgewidth=0, color='black', markersize=3)
            ax[i].plot(n[ok == False].flatten(),noise[ok == False].flatten(), marker='o', linewidth=0, alpha=0.02, color='red', markersize=3, markeredgewidth=0)
            bx, by, be= zachopy.oned.binto(x=n, y=noise, binwidth=1000, yuncertainty=None, robust=True, sem=True)
            ax[i].errorbar(bx, by, be, color='orange', linewidth=5, elinewidth=5, capthick=5, alpha=0.5)

            x = np.linspace(0, n[ok].max(), 100)
            ax[i].plot(x, noisemodel(x,gain,readnoise,lamp), color='green', linewidth=5, alpha=1.0)
            ax[i].plot(x, noisemodel(x,1e20,readnoise,0.0), color='green', linestyle='--',linewidth=5, alpha=0.2)
            ax[i].plot(x, noisemodel(x,gain,0.0, 0.0), color='green', linestyle='--',linewidth=5, alpha=0.2)
            ax[i].plot(x, noisemodel(x,1e20, 0.0, lamp), color='green', linestyle='--',linewidth=5, alpha=0.2)

            ax[i].set_title('c{0} \n gain = {1:.2f} \n readnoise = {2:.2f} \n lamp variability = {3:.3f}'.format(i+1, gain, readnoise, lamp))
            ax[i].set_xlim(np.min(n[ok]), np.max(n[ok]))
            ax[i].set_ylim(np.maximum(np.min(noise[ok]), 1.0), np.max(noise[ok]))
            ax[i].set_yscale('log')
            ax[i].set_xscale('log')
            plt.draw()
            #ax[i].plot(n[ok].flatten(), deviates(fit.params, n=n[ok], noise=noise[ok])[-1].flatten(), marker='o', linewidth=0, alpha=0.005, color='black', markersize=3)
            #ax[i].plot(n[ok==False].flatten(), deviates(fit.params, n=n[ok == False], noise=noise[ok == False])[-1].flatten(), marker='o', linewidth=0, alpha=0.1, color='red', markersize=3)
            assert('n' not in self.input('Does this gain estimate seem reasonable? [Y,n]').lower())

          np.savetxt(self.obs.workingDirectory + 'gains.txt', gains)
          self.gains = gains

  def createStackedImage(self, n, visualize=True, imageType=None, threshold=5.0, truncation=100):
    '''Take an outlier-rejected stack of a series of images (requires enough memory to hold them all.)'''

    stride = np.maximum(len(n)/truncation, 1)

    if stride > 1:
        self.speak('stacking {0}/{2} {1} images'.format(len(n),imageType,truncation))
    else:
        self.speak('stacking {0} {1} images'.format(len(n), imageType))

    array = self.ccd.loadImages(n[::stride], imageType=imageType)
    if len(array.shape) <=2:
      return array, array*0


    median, noise = zachopy.twod.stack(array,axis=0,threshold=threshold)

    #self.ccd.display.many(array, depth=0, clobber=True)
    return median, noise

  def createMeanImage(self, n, cosmic=True, visualize=False, imageType=None):
    '''Take the mean of a series of images, less memory intensive than median.'''
    stride = np.maximum(len(n)/100, 1)
    for i in np.arange(0, len(n), stride):
      data = self.ccd.readData(n[i], imageType=imageType)
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
    badPixelFilename = self.obs.workingDirectory + 'master_BadPixels.fits'
    try:
      self.images['BadPixels'] = readFitsData(badPixelFilename)
      self.speak( "loaded bad pixel mask from {0}".format(badPixelFilename))
    except:
      self.speak( "creating bad pixel mask from the master flat frames")
      c = self.ccd#CCD(self.obs, calib=self)

      cube = []
      for n in self.obs.nWideFlat:
        c.set(n, 'WideFlat')
        cube.append(c.readData())

      cube = np.array(cube)

      median = np.median(cube,0)
      writeFitsData(median, 'median.fits')
      noise = np.median(np.abs(cube - median.reshape(1,cube.shape[1], cube.shape[2])), 0)
      plt.figure('bad pixel mask')
      ax = plt.subplot()
      ax.plot(median.flatten(), noise.flatten(), color='black', alpha=0.5, marker='o', markersize=4, markeredgewidth=0, linewidth=0)
      ax.set_yscale('log')
      if self.obs.instrument == 'LDSS3C': bad = (noise < 0.05*np.sqrt(median)) | (noise == 0) | (median == 0) | (median < 0) | (self.bias() > 10000) | (self.dark() > 100)
      elif self.obs.instrument == 'IMACS': bad = (noise < 0.05*np.sqrt(median)) | (noise == 0) | (median == 0) | (median < 0)
      ax.plot(median[bad].flatten(), noise[bad].flatten(), color='red', alpha=0.5, marker='o', markersize=10, markeredgecolor='red', linewidth=0)
      ax.set_xlabel('Fluence')
      ax.set_ylabel('RMS')
      self.images['BadPixels'] = bad.astype(np.int)
      if visualize:
        self.display.one(self.images['BadPixels'])
        answer = self.input("Does the bad pixel mask seem reasonable? [Y,n]").lower()
        assert('n' not in answer)
      writeFitsData(self.images['BadPixels'], badPixelFilename)




  def createMasterImage(self, imageType=None, remake=False):

    print imageType

    # if no name included, do nothing
    if imageType is None:
      return

    # set the CCD to a particular image type
    self.ccd.set(n=None, imageType=imageType)

    # make sure this is a master image that is possible to make
    assert(imageType in self.obs.cal_dictionary.keys())

    noisestring = 'StdDev'
    self.speak( "populating the master {0} image".format(imageType))
    masterFilePrefix = self.obs.workingDirectory + "master_{0}".format(imageType)
    try:
        self.images[imageType] = readFitsData(masterFilePrefix + '.fits')
        self.images[imageType+noisestring] = readFitsData(masterFilePrefix + noisestring + '.fits')
        self.speak( "loaded {0} from {1}.fits".format(imageType, masterFilePrefix))
    except IOError:
        self.speak("creating from images " + zachopy.utils.truncate( str( self.obs.cal_dictionary[imageType]), n=30))
        self.images[imageType], self.images[imageType+noisestring] = self.createStackedImage(self.obs.cal_dictionary[imageType], imageType=imageType)
        writeFitsData(self.images[imageType], masterFilePrefix  + '.fits')
        writeFitsData(self.images[imageType+noisestring],masterFilePrefix + noisestring + '.fits')
        #self.display.one(self.images[imageType+noisestring], clobber=True)
        #self.display.one(self.images[imageType], clobber=False)
        #self.display.single()
        #self.display.zoom()
        #self.display.scale('log', limits=[0,np.percentile(self.images[imageType],99)])
        #assert('n' not in self.input("Do you like master image {0}? [Y,n]".format(imageType)).lower())

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
        if 'Science' not in k:
            self.createMasterImage(k, remake=remake)
