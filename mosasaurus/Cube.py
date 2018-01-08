from .imports import *
from .Observation import Observation
from numpy.polynomial import Legendre

cubehelp = '''
Independental measured variables:

=======================
stored in cube.temporal
=======================
airmass(time)
this is the overall average airmass of the field

rotatore(time)
this the instrument rotator angle -- it probably tells us something about the
instrument's changing illumination and flexure.

=======================
stored in cube.squares
=======================

centroid(star,time)
the centroid of the star in the cross-dispersion direction. this variable is the
median of the centroids across all wavelengths, representing an overall shift of
the star's position through the spectrograph's optics

width(star,time)
the width of the star in the cross-dispersion direction. this variable is the
median of the widths across all wavelengths, representing an overall change,
most likely due to seeing.

shift(star,time)
by how much did we have to shift the star in the wavelength direction, to
make its spectral features line up with a reference spectrum. as the zeropoint
of the wavelength calibration is set by the position of the star in the focal
plane, this tracks the motion of the star in the dispersion direction (and
probably other stuff too)

=======================
stored in cube.cubes
=======================

raw_counts(star,time,wavelength)
the flux from each star, integrated over the cross-dispersion direction and with
flux from the diffuse sky subtracted. this is what we use to make light curves!

sky(star,time,wavelength)
the extracted sky brightness, which was subtracted in an average sense, from the
the 1D stellar spectrum. if this value was estimated incorrectly, it could
influence the chromatic light curves.

dcentroid(star,time,wavelength)
the centroid of the star in the cross-dispersion direction, measured at a
particular wavelength, relative to "centroid" (see above). this represents
changes in the position on the detection caused by things like internal flexure
or atmospheric refraction.

dwidth(star,time,wavelength)
the width of the star in the cross-dispersion direction, measured at a
particular wavelength, relative to "width". flux correlations with dwidth
could be diagnostic of problems with the sky subtraction, or the extraction

peak(star,time,wavelength)
the brightness of the brightest pixel in the cross-dispersion direction,
at each wavelength, for each star. it will correlate strongly with seeing,
but is a slightly different tracer. correlations with "peak" that cannot be
explained by seeing alone might point to problems in the detector non-linearity
'''


plt.ion()

# set a cmap for stars
starcm = zachopy.cmaps.one2another('magenta', 'limegreen')

class Cube(Talker):
  '''
  Cube object stores a time-wavelength-flux datacube.

  It has one unique extraction width associated with it,
  even though the individual supersampled spectra from
  which it's been compiled may have many widths available.
  '''

  # keep - initialization
  def __init__(self, obs, width=None, remake=False, max=None, **kwargs):
    '''Initialize a data cube and populate it with data.

    You need to supply an Observation object, and one extraction width.
    (NOTE, this will break when using different extraction widths for
    different stars within the same observation).

    '''
    Talker.__init__(self, line=200, **kwargs)


    # set up the basics
    self.obs = obs
    self.reducer = self.obs.reducer

    # the temporary filename for the cube
    self.tempfilename = os.path.join(self.directory, 'tempSpectralCube.npy')
    self.cubekeys = ['raw_counts', 'sky',  'centroid', 'width', 'peak']

    # what attributes are savable/loadable?
    self.savable = ['cubes', 'squares', 'temporal', 'spectral', 'stellar']

    # define the dimensions of a cube
    self.sindex = 0
    self.tindex = 1
    self.windex = 2

    # set initial values for these (they may change)
    self.numberofstars = len(self.starDirectories)
    self.numberoftimes = len(self.obs.fileprefixes['science'])

    # make sure the user has specified a width
    if width == None:
        self.speak('Yikes! Nobody set an extraction width with which to create this cube!')
        self.listWidths()
        width = np.float(self.input("What width would you like? [as a number]"))

    self.width = width

  def __repr__(self):
      return '<Cube of {} Spectra for {} Stars for {}>'.format(self.numberoftimes, self.numberofstars, self.obs)

  # keep - shortcut for this cube's directory
  @property
  def directory(self):
    return self.reducer.extractionDirectory

  # keep - lists the stars that are available (and shows finder chart)
  def listStars(self):
    '''
    List the available stars for this cube,
    and open the finder chart.
    '''
    self.speak('the stars available for {} are'.format(self.obs))
    for s in self.stars:
        self.speak("   '{}'".format(s))
    self.speak('as shown in the finder chart')
    os.system('open {}'.format(os.path.join(self.directory, 'genericfinderchart.pdf')))

  # keep - setting which stars to use as comparisions
  def setStars(self, target=None, comparisons=[1], other=None):
    '''
    Decide with stars should be used as
    the target, the comparison[s], or something else.
    '''

    if target == None:
        self.listStars()
        raise ValueError("Please specify the target and comparison(s).")

    # this should be just one star (with a string name)
    self.target = target

    # these can be one or more stars (and we'll make sure they're 1D arrays)
    self.comparisons = np.array([comparisons]).flatten()
    self.other = np.array([other]).flatten()

    self.speak('setting the [target] star to {}'.format(self.target))
    self.speak('setting the [comparisons] star(s) to {}'.format(self.comparisons))
    self.speak('setting the [other] stars to {}'.format(self.other))

  # keep - helps figure out options for widths to use
  def listWidths(self):
    '''
    Print out the widths associated with each star.
    '''
    self.speak('checking widths available for {}'.format(self.obs.reducer.extractionDirectory))
    for s in self.starDirectories:
      extractionandaperture = s.split('extraction_')[-1]
      self.speak(' [{}] contains:'.format(extractionandaperture))
      files = glob.glob(os.path.join(s, 'supersampled_*.npy'))
      supersampled = np.load(files[0])[()]
      for k in supersampled.keys():
        if 'raw_counts' in k:
            width = k.split('_')[-1]
            self.speak('   {}'.format(width))

  # keep - allow us to populate a cube
  def populate(self, remake=False, max=None, visualize=True, shift=True):
    '''
    Populate this cube with data, either by loading a saved cube,
    or by loading all of the spectra individually into a new cube.
    '''

    self.shift=shift
    try:
        # does the cube already exist?
        self.cubes
        # does its width already match the desired width?
        assert(self.width == self.meta['width'])
        self.speak('the cube for {} has already been loaded'.format(os.path.split(self.filename)[-1]))
    except (AttributeError,AssertionError):
        try:
            # we might be forced by a keyword to remake
            assert(remake==False)
            # if a single cube file already exists, let's load that
            self.load()
        except (IOError, AssertionError): #AttributeError
            # otherwise, load individual spectra an populate this cube
            self.loadSpectra(remake=remake, max=max, visualize=visualize)

    # make some summary images of this cube
    #self.imageCube()

  # keep - shortcut for the star directories
  @property
  def starDirectories(self):
      '''return the absolute paths of all star directories created for this mask'''
      return glob.glob(os.path.join(self.directory, 'aperture_*'))

  # keep - the filename of this cube, based on lots of parameters
  @property
  def filename(self):
    '''the filename of this cube (depends on whether it is shifted, or not)'''
    s = {True:'shifted', False:'raw'}
    name = 'spectralCube_{}_{}_{}stars_{}spectra_{}_{}.npy'.format(
        self.obs.target.name,
        self.obs.night.name,
        self.numberofstars,
        self.numberoftimes,
        self.widthkey.replace('.0', ''),
        s[self.shift])
    return os.path.join(self.directory, name)

  # keep - a short cut for the stars that are present in the cube
  @property
  def stars(self):
      '''return a list of star names that belong in this cube'''
      return self.stellar['aperture']

  # keep - a short cut for the width key (associated with a numerical extraction width)
  @property
  def widthkey(self):
    return '{:04.1f}px'.format(self.width)

  # keep - load individual spectra, and populate the cube!
  def loadSpectra(self, remake=False, visualize=True, max=None):
    """
    The opens the supersample spectra in the extraction directory,
    one by one, and uses them to populate the entries in the cube,
    for the currently specified extraction width.
    """

    # 3d, stars x time x wavelength
    self.cubes = {}
    # 2d, stars x time
    self.squares = {}
    # 1d, time
    self.temporal = {}
    # 1d, wavelength
    self.spectral = {}
    # 1d, star
    self.stellar = {}
    # other details about this extraction
    self.meta = {}

    # update
    self.speak("Loading the spectral cube.")

    # define the number of stars and times we're looking for
    self.numberofstars = len(self.starDirectories)
    self.numberoftimes = len(self.obs.fileprefixes['science'])
    if max is not None:
        self.numberoftimes = max
    truncate = False

    # load the headers (from the observation object)
    self.headers = self.obs.headers

    # load the names of the stars
    self.stellar['aperture'] = [x.split('/')[-1] for x in self.starDirectories]

    # loop over the spectra
    for timepoint in range(self.numberoftimes):
        # pull out the file prefix for this star
        fileprefix = self.obs.fileprefixes['science'][timepoint]

        # loop over all the stars
        for istar, star in enumerate(self.stars):

          # find the available spectrum

          # use the updated wavelengths if we're looking at the shifted cube
          if self.shift:
            extractedFile = os.path.join(self.starDirectories[istar], 'updatedwavelengthextracted_{0}.npy'.format(fileprefix))
            spectrumFile = os.path.join(self.starDirectories[istar], 'updatedwavelengthsupersampled_{0}.npy'.format(fileprefix))
          else:
            extractedFile = os.path.join(self.starDirectories[istar], 'extracted_{0}.npy'.format(fileprefix))
            spectrumFile = os.path.join(self.starDirectories[istar], 'supersampled_{0}.npy'.format(fileprefix))

          self.speak('trying to load {0}'.format(spectrumFile))
          # load the extracted spectrum (or truncate the cubes at this point)
          try:
              supersampled = np.load(spectrumFile)[()]
              self.speak('loaded {0}'.format(spectrumFile))
              extracted = np.load(extractedFile)[()]
              self.speak('loaded {0}'.format(extractedFile))
          except IOError:
              # if we've run out of spectra to load, then truncate
              truncate = True
              self.speak('failed to find {}'.format(spectrumFile))
              self.speak('truncating cube!')
              break

          try:
            # have I already loaded these ingredients?
            self.spectral['wavelength']
            self.spectral['fractionofapixel']
            self.numberofwavelengths
          except (KeyError,AttributeError):
            # define some useful arrays
            self.spectral['wavelength'] = supersampled['wavelength']
            self.spectral['fractionofapixel'] = supersampled['fractionofapixel']
            self.numberofwavelengths = len(self.spectral['wavelength'])

          # make sure the wavelength grid matches what we've stored (should be same across all stars)
          assert((self.spectral['wavelength'] == supersampled['wavelength']).all())

          # loop over the measurement types and populate the cubes
          for key in self.cubekeys + ['ok']:

            # make sure a cube exists for this key
            try:
                self.cubes[key]
            except KeyError:
                self.cubes[key] = {}

            # make sure a cube entry exists for this star (an array of times and wavelengths)
            try:
              self.cubes[key][star]
            except KeyError:
              if key == 'ok':
                self.cubes[key][star] = np.ones((self.numberoftimes, self.numberofwavelengths)).astype(np.bool)
              else:
                self.cubes[key][star] = np.zeros((self.numberoftimes, self.numberofwavelengths)).astype(np.float32)
            self.speak("updating cubes['{key}']['{star}'][{timepoint},:]".format(**locals()))

            # populate with the supersampled spectrum
            if key != 'ok':
                self.cubes[key][star][timepoint,:] = supersampled[key + '_' + self.widthkey]

            if 'raw_counts' in key:
                s = sum(self.cubes[key][star][timepoint,:])
                self.speak('(raw_counts sum to {} for {})'.format(s, fileprefix))
                assert(s>0.0)

          # pull out data from the (unsupersampled) spectra to populate a square with dimensions self.numberofstars x self.numberoftimes
          for key in ['sky', 'width', 'centroid']:#, 'cosmicdiagnostic']:

              try:
                  self.squares[key]
              except KeyError:
                  self.squares[key] = {}
              try:
                  self.squares[key][star]
              except KeyError:
                  self.squares[key][star] = np.zeros(self.numberoftimes).astype(np.float32)

              self.squares[key][star][timepoint] = np.nanmedian(extracted[self.width][key])
              self.speak("updating squares['{key}']['{star}'][{timepoint}]".format(**locals()))

        # if we've run out of spectra, then break out of the loop (with truncated cubes)
        if truncate:
            break

        self.speak('{0}/{1} spectra loaded into cube'.format(timepoint, self.numberoftimes))

        # if the spectra for all stars were successfully loaded, then
        try:
            self.temporal['fileprefix']
        except KeyError:
            self.temporal['fileprefix'] = []
        self.temporal['fileprefix'].append(fileprefix)

    # make sure everything is truncated properly
    if truncate:
        self.speak("couldn't find all requested spectra, so truncated cube at a length of {0}".format(timepoint))
        for key in self.cubes.keys():
            self.cubes[key] = self.cubes[key][star][0:timepoint,:]
        for key in self.squares.keys():
            self.squares[key] = self.squares[key][star][0:timepoint]

    # keep track of purely time-dependent quantities
    self.temporal = astropy.table.Table(self.headers)[0:self.numberoftimes]
    self.temporal['ok'] = np.ones(self.numberoftimes).astype(np.bool)#self.temporal['cosmicdiagnostic'] < self.obs.cosmicAbandon

    # store some metadata
    self.meta['width'] = self.width
    self.meta['target'] = self.obs.target.name
    self.meta['night'] = self.obs.night.name
    self.meta['instrument'] = self.obs.instrument.name
    self.meta['extractiondefaults'] = self.obs.instrument.extractiondefaults


    if self.shift:
        self.shiftCube(plot=visualize)

    self.speak("Done loading spectral cube.")

    #self.markBad()
    self.save()

  # keep - save the entire cube to a single file
  def save(self):
      self.speak('attempting to save the cube of loaded, shifted, compiled spectra to {0}'.format(self.filename))
      tosave = {}
      for thing in self.savable:
          tosave[thing] = self.__dict__[thing]
          self.speak('  including [{0}] in the saved cube structure'.format(thing))
      np.save(self.filename, tosave)

  # keep - load the entire cube from a single file
  def load(self):
      self.speak('attempting to load previously saved cubes from...')
      self.speak('{0}'.format(self.filename))
      loaded = np.load(self.filename)[()]
      for thing in self.savable:
          self.__dict__[thing] = loaded[thing]
          self.speak('  loading [{0}] from the saved cube structure'.format(thing))
      self.numberofstars = len(self.stellar['aperture'])
      self.numberofwavelengths = len(self.spectral['wavelength'])
      self.numberoftimes = len(self.temporal)

  # keep - automated tool to package this observation into a nice directory
  def packageandexport(self):
      '''make a tidy package that contains necessary information for sending to someone else'''

      # create a directory for this packaged cube to land
      base = self.obs.reducer.extractionDirectory
      toexport = os.path.join(base, 'cube_' + self.meta['target'] + '_' + self.meta['night'])
      mkdir(toexport)
      commandstorun = []

      # copy the cube npy file
      commandstorun.append('cp {} {}/.'.format(self.filename, toexport))

      # copy the finder chart
      commandstorun.append('cp {} {}/.'.format(os.path.join(base, 'genericfinderchart.pdf'), toexport))

      # copy the finder chart
      commandstorun.append('cp {} {}/.'.format(os.path.join(base,'extractionCenters.txt'), toexport))

      # loop over the star directories
      for s in self.stars:
        stardir = os.path.join(toexport, s)
        mkdir(stardir)
        pdfs = glob.glob(os.path.join(base, s, '*aperture*.pdf'))
        for p in pdfs:
            # copy the finder chart
            commandstorun.append('cp {} {}/.'.format(p, stardir))

        apertureimage = os.path.join(base, s, 'animatedextraction/formovie_00000.png')
        commandstorun.append('cp {} {}/{}_whichaperture.png'.format(apertureimage, stardir, s))

      for c in commandstorun:
          print(c)
          os.system(c)

  @property
  def cubeMegaComparison(self):
    '''
    Return a fake cube entry for an combined mega-comparison star.
    '''
    if len(self.comparisons) == 1:
        d = {}
        star = self.comparisons[0]
        for key in self.cubes.keys():
            d[key] = self.cubes[key][star]
        return d
    else:
        raise NameError("Darn it -- the mega-comparison hasn't been implemented yet!")

  # keep - makes tidy image plots of spectral quantities
  def imageCube(self, normalized=False, remake=False):
    '''
    Show an imshow of every cube key.
    '''

    # we'll plot various types of images
    options = {
                'raw':'Raw Extracted Quantities',
                'wavelength':'Wavelength-normalized Quantities',
                'comparison':'Comparison-divided Quantities'
              }

    for option, description in options.items():
        filename = os.path.join(self.directory, 'imagedcube_{}.pdf'.format(option))
        if os.path.exists(filename) and (remake == False):
            continue
        nrows = len(self.cubekeys)
        ncols = len(self.stars)

        fi, ax = plt.subplots(nrows, ncols,
                                sharex=True, sharey=True, figsize=(12,8),
                                gridspec_kw=dict(hspace=0.1, wspace=0.02))
        plt.suptitle('{}, [width={}]\n{}'.format(description, self.widthkey, self.obs))

        w = self.spectral['wavelength']
        t = np.arange(self.numberoftimes)

        # set up the imshow parameters
        imkw = dict(
            extent = [np.min(w), np.max(w), np.min(t), np.max(t)],
            cmap = 'gray',
            interpolation='nearest',
            aspect = 'auto',
            origin = 'lower'
        )

        for i, key in enumerate(self.cubekeys):
            for j, star in enumerate(self.stars):
                a = ax[i,j]


                if option == 'raw':
                    # don't modify the measurements at all
                    self.speak('displaying the raw measurements for {}'.format(key))
                    z = self.cubes[key][star]

                if option == 'wavelength':
                    # normalize along the wavelength axis
                    self.speak('displaying the measurements normalized by their median spectrum for {}'.format(key))
                    z = self.cubes[key][star]
                    oned = np.median(z, 0)
                    z = z/oned[np.newaxis,:]

                if option == 'comparison':
                    # divide by the comparison star[s]

                    self.speak('displaying the measurements divided by the comparison[s] for {}'.format(key))
                    target = self.cubes[key][star]
                    comparison = self.cubeMegaComparison[key]
                    z = target/comparison
                    oned = np.median(z, 0)
                    z = z/oned[np.newaxis,:]

                vmin, vmax = np.percentile(z, [1,99])
                a.imshow(z, vmin=vmin, vmax=vmax, **imkw)

                # fuss with the axis labels
                if j == 0:
                    a.set_ylabel('{}\n(timepoints)'.format(key))
                else:
                    plt.setp(a.get_yticklabels(), visible=False)
                if i == 0:
                    a.set_title('{}'.format(star))

                if i == (len(self.cubekeys)-1):
                    a.set_xlabel('Wavelength (angstroms)')
                else:
                    plt.setp(a.get_xticklabels(), visible=False)

        filename = os.path.join(self.directory, 'imagedcube_{}.pdf'.format(option))
        plt.savefig(filename)
        self.speak('saved image of this cube to {}'.format(filename))


  def nudgeWavelengths(self):
    '''
    This function loops through the extracted*.npy files
    and determines corrected wavelength solutions that will match up
    their lines. It writes a new array of wavelengths for each,
    into a file that aligns with the original extracted file.

    (Can be run before an unshifted cube is populated.)

    The code for this was developed by Hannah Diamond-Lowe.

    '''

    def align_lines(stars, line_range, offset, line_pos, line_poses, plot=False):
        '''
        This function takes
            stars:    the extracted*.npy file
            line_range:    the pre-determined fixed alignment range for a given feature to use for the alignment
            offset:   how many angstroms to shift over when performing the cross-correlation of a given feature in a spectrum
            line_pos:    a list that this function will append values to; amount of shift in wavelength space for each star width
            line_poses:    a list of line_pos lists; amount of shifts in wavelength space for each star with all its widths
            plot:  do you want plots to pop up? (T/F)

        This code assumes that the stellar spectra may shift slowly throughout the night but do not suddenly jump around from one exposure to the next.

        '''

        #### ZKBT: pulling out starmaster from the end of stars?
        starmaster = stars[-1]

        # find the wavelengths between the wavelength range of the feature for the master star
        idxmaster = (starmaster['wavelength']>=line_range[0])*(starmaster['wavelength']<=line_range[1])
        # get the corresponding values from raw_counts; this will be the main thing to correlate against
        corrmaster = starmaster[width]['raw_counts'][idxmaster]

        ####### KLUDGE?????
        # corrmaster = zachopy.oned.subtractContinuum(corrmaster, n=2) + 0.0

        # this is where the "true" value of the line is; i.e., the reference position on the master spectrum
        line = starmaster['wavelength'][idxmaster][0]
        # where that line falls in pixel space; uses the wavelength solution from mosasaurus
        linepx = np.interp(line, starmaster['wavelength'], starmaster['w'])

        # list of stars is the stars from "aperture" in mosasaurus, plus the master star appended at the end
        for s in range(len(stars)):

            if (plot) and (s != (len(stars)-1)): self.speak('checking correlations for {}'.format(self.stars[s]))

            # find the wavelengths between the wavelenth range of the feature for this particular star
            idxstar = (stars[s]['wavelength']>=line_range[0])*(stars[s]['wavelength']<=line_range[1])
            # this is the spectrum of the star that you want to correlate to the master
            corrstar = stars[s][width]['raw_counts']
            ####### KLUDGE?????
            # corrstar = zachopy.oned.subtractContinuum(corrstar) + 0.0

            # need to know now many wavelengths are being covered
            arraylen = len(np.where(idxmaster)[0])
            # this is the pixel where we will start the shift
            initpx = np.where(idxstar)[0][0] - offset
            #if corrstar[initpx] == 0.0:
            #    initpx = np.where(corrstar != 0.)[0][0]
            #    corrmaster = wavemaster[initpx+offset:initpx+offset+arraylen]
            # empty lists for the correlation coeffients and the actual pixel shifts
            corrs, shifts = [], []

            # set the number of shifts you will try in the cross correlation; the finer the shift the more accurate everything will be but it will take longer
            for shift in np.linspace(0, arraylen+offset, (arraylen+offset)*10+1):
                # this is to make sure we're not losing counts at the ends of pixels when we're shifting by sub-pixels
                newbit = shift%1
                therest = 1 - newbit
                startpx = int(np.floor(shift))
                # create the array at the sub-pixel level that will be compared to the master star feature
                corrarray = corrstar[initpx+startpx : initpx+startpx+arraylen]*therest + corrstar[initpx+startpx+1 : initpx+startpx+arraylen+1]*newbit
                #if shift == 0:
                #    plt.plot(corrmaster)
                #    plt.plot(corrarray)
                #    plt.show()
                corrs.append(np.corrcoef(corrmaster, corrarray)[0,1])
                shifts.append(corrarray)

            if plot == True and (s != len(stars)-1):
                # can inspect where the code thinks the peak of the correlation is (careful! it doesn't always pick the right one!)
                #print('1st corr')
                color = self.starcolor(self.stars[s])
                plt.plot(corrs, color=color)
                plt.axvline(np.where(corrs == np.max(corrs))[0], color=color)
                plt.title(str(line_range))
                plt.xlim(0,400)
                #plt.draw()
                #plt.show()
                #a = input('hmmm?')

            # try a shift based on the correlation coefficients
            # this first try may be wrong so we have to compare to past shifts
            firstshiftind = np.where(corrs == np.max(corrs))[0][0]
            firstrange = np.linspace(0, arraylen+offset, (arraylen+offset)*10+1)
            firstshift = firstrange[firstshiftind]
            # offset should be reasonably small since we assume the spectral features do not suddenly jump around from one exposure to the next
            if len(line_poses) < 13:
                while ((firstshift - offset) > 5.) or ((firstshift - offset) < -5.):
                    self.speak('doing an extra correlation {}: {}'.format(line_range, firstshift-offset))
                    currentmaxind = np.where(corrs == np.max(corrs))
                    # delete the correlation coefficient that is the maximum but is not providing a reasonably small shift
                    corrs = np.delete(corrs, currentmaxind)
                    # find a new shift value
                    firstshiftind = np.where(corrs == np.max(corrs))[0][0]
                    firstshift = firstrange[firstshiftind]
                    if ((firstshift - offset) < 5.) and ((firstshift - offset) > -5.) and plot == True and s != len(stars)-1:
                        color = self.starcolor(self.stars[s])
                        plt.plot(corrs, color=color)
                        plt.axvline(np.where(corrs == np.max(corrs))[0], color=color)
                        assert(False)
                        #plt.title(str(line_range))
                        #plt.show()

            # once we have made enough reasonably small shifts we can use past shifts to determine whether or not the next shift is a reasonable jump
            # there is probably a better way to do this... like going through the whole thing and then looking for outliers
            else:
                # look at the last 13 shifts in line position for this particular star
                med = np.median(np.array(line_poses)[:, s][-13:] - np.array(line_poses)[:, -1][-13:])
                std = np.std(np.array(line_poses)[:, s][-13:] - np.array(line_poses)[:, -1][-13:])
                # this is apparently a good number of sigma to go by
                timesstd = 22
                # eventually the value that gets save in line_poses is a wavelength; make the transformation into wavelength space so that we can compare the current shift to the previous ones
                if (pxtowavemaster(linepx + firstshift - offset)-line > (med + timesstd*std)) or (pxtowavemaster(linepx + firstshift - offset)-line < (med - timesstd*std)):
                    print('range: ', med - timesstd*std, '-', med + timesstd*std, pxtowavemaster((linepx + firstshift - offset))-line)
                    print('stddev corr', line_range)
                    corrind = np.where((pxtowavemaster(linepx + firstrange - offset)-line >= (med - timesstd*std)) & (pxtowavemaster(linepx + firstrange - offset)-line <= (med + timesstd*std)))
                    corrsclipped = np.array(corrs)[corrind]
                    firstshiftind = np.where(corrsclipped == np.max(corrsclipped))[0][0] + corrind[0][0]
                    firstshift = firstrange[firstshiftind]
                    #if plot == True and s != len(stars)-1:
                    #    #plt.figure('something - ' + str(line_range))
                    #    plt.plot(corrs, color='orange')
                    #    plt.axvline(np.where(corrsclipped == np.max(corrsclipped))[0][0] + corrind[0][0], color='orange')

            # interpolate shift using a parabola to find the true peak cross correlation value
            parabinds = range(firstshiftind-2,firstshiftind+3)
            parabrange = np.linspace(parabinds[0], parabinds[-1], 100*(len(parabinds)-1))
            parabfit = np.polyfit(parabinds, np.array(corrs)[parabinds], 2)
            parabfit1d = np.poly1d(parabfit)
            parabcorrs = parabfit1d(parabrange)
            parabshiftind = np.where(parabcorrs == np.max(parabcorrs))[0]
            parabshift = parabrange[parabshiftind]
            #if plot == True:
            #    plt.plot(parabinds, np.array(corrs)[parabinds])
            #    plt.plot(parabrange, parabcorrs)
            #    plt.axvline(parabshift)
            #   plt.show()
            frac = parabshift - firstshiftind
            finalshift = firstshift + frac*np.median(np.diff(firstrange))

            # transform all these shifts in pixel space into wavelength space
            newlinepx = (linepx + finalshift - offset)[0]
            newline = pxtowavemaster(newlinepx)
            #print newline
            # remember that the last "star" in the list is actually just the master spectrum; again, not sure why I did it this way
            if s == len(stars)-1: line_pos.append(line)
            else: line_pos.append(newline)
            if plot == True: self.speak(' shift = {}'.format(finalshift - offset))
            #print newline

            # this was left-over from some test I did to see how much the GJ1132 spectrum was stretching/shifting throughout the night
            #shift1132.append((starline+finalshift-offset)[0])

            #newbit = finalshift%1
            #therest = 1 - newbit
            #startpx = int(np.floor(finalshift))
            #newarray = corrstar[initpx+startpx:initpx+startpx+arraylen]*therest + corrstar[initpx+startpx+1:initpx+startpx+arraylen+1]*newbit

            #if plot == True:
            #    plt.plot(newarray, label=str(s))

        #if plot == True:
        #    plt.plot(corrmaster, label='corrmaster')
        #    plt.title(str(line_range))
        #    plt.legend()
        #    plt.show()

    # pick one star and exposures that will be the reference for all
    # what's the master's directory
    masterstar = self.stars[0]
    masterexposure = self.obs.fileprefixes['science'][0]
    aperturedirectory = os.path.join(self.directory, masterstar)
    # what's the master star file?
    starmasterstr = os.path.join(aperturedirectory,'extracted_{}.npy'.format(masterexposure))
    starmaster = np.load(starmasterstr)[()]
    self.speak('loaded {} as the master reference'.format(masterstar))

    # load the original wavelength calibration
    wavecalfile = os.path.join(aperturedirectory, '{}_wavelengthcalibration.npy'.format(masterstar))
    coef, domain =  np.load(wavecalfile)[()]
    self.speak('loaded {} as the original wavelength calibration file for {}'.format(wavecalfile, masterstar))


    # reacreate the wavelength solution (way to go from pixel space to wavelength space) from mosasaurus
    pxtowavemaster = Legendre(coef, domain)
    # apertures form a given night (Need to have run mosasaurus in ipython or similar and then copy and paste this monster code in. This is not ideal.)
    makeplot = False
    ### UV_poses, O2_poses, Ca1_poses, Ca2_poses, Ca3_poses, H2O_poses = [], [], [], [], [], []
    x_poses = []
    #shifts_1132 = []

    # the alignment ranges are set by the instrument and gratings
    alignmentranges = self.obs.instrument.alignmentranges


    plt.clf()

    rangenames = alignmentranges.keys()
    all_poses = {}
    for k in alignmentranges.keys():
        all_poses[k] = []
    # Fixed alignment ranges for each prominent freature
    ### align_UV = (6870, 6900)
    ### align_O2 = (7580, 7650)
    ### #align_H2Osmall = (8220, 8260)
    ### align_Ca1 = (8490, 8525)
    ### align_Ca2 = (8535, 8580)
    ### align_Ca3 = (8650, 8700)
    ### align_H2O = (9300, 9700)
    for prefix in self.obs.fileprefixes['science']:
        self.speak('checking line positions in {}'.format(prefix))
        stars = []
        for d in self.starDirectories:
            extractedpathname = os.path.join(d, 'extracted_{}.npy'.format(prefix))
            stars.append(np.load(extractedpathname)[()])
        # append to the list of stars the master spectrum; not sure why I did it this way but we do need to know the "master" line position in wavelength space so we can later calculate how many wavelengths to shift over
        stars.append(starmaster)

        # loop over all the widths
        widths = [x for x in stars[0].keys() if type(x) != str]
        for width in widths:
            # make a directory
            corrdir = os.path.join(self.directory, 'correlations_{}'.format(width))
            mkdir(corrdir)

            self.speak('  aligning lines for [{:04.1f}px]'.format(width))

            # the change in position of each of these features
            ### UV_pos, O2_pos, Ca1_pos, Ca2_pos, Ca3_pos, H2O_pos = [], [], [], [], [], []
            # NOW GO UP TO THE aling_star FUNCTION
            #align_lines(stars, align_UV, 10, UV_pos, UV_poses, makeplot)
            #align_lines(stars, align_O2, 20, O2_pos, O2_poses, makeplot)

            # a dictionary to store the pos for each line
            all_pos = {}
            for k in rangenames:
                all_pos[k] = []

            if True:
                plt.clf()
                plt.figure('correlations', figsize=(20,3), dpi=30)
                gs = plt.matplotlib.gridspec.GridSpec(1, len(rangenames),
                            hspace=0.2, wspace=0.25, left=0.05, right=0.95)
                self.ax_corr = {}

            for i, k in enumerate(rangenames):

                self.ax_corr[k] = plt.subplot(gs[i])
                plt.xlabel('Offset')

                temp = all_pos[k]
                align_lines(stars, alignmentranges[k], 10, temp, all_poses[k], True)
                all_pos[k] = temp
                all_poses[k].append(all_pos[k])


            pltfilename = os.path.join(corrdir, 'wavelengthcorrelations_{}.pdf'.format(prefix))
            plt.savefig(pltfilename)
            self.speak('saved plot to {}'.format(pltfilename))

            ### align_lines(stars, align_Ca1, 10, Ca1_pos, Ca1_poses, makeplot)
            ### align_lines(stars, align_Ca2, 10, Ca2_pos, Ca2_poses, makeplot)
            ### align_lines(stars, align_Ca3, 10, Ca3_pos, Ca3_poses, makeplot)
            ### align_lines(stars, align_H2O, 10, H2O_pos, H2O_poses, makeplot)

            #UV_poses.append(UV_pos)
            ### O2_poses.append(O2_pos)
            ### Ca1_poses.append(Ca1_pos)
            ### Ca2_poses.append(Ca2_pos)
            ### Ca3_poses.append(Ca3_pos)
            ### H2O_poses.append(H2O_pos)


            # do some list re-arranging; we need to know what the shift and stretch in wavelength are for each star
            # x_pos just means any of the feature positions we've calculated
            # x_pos = np.array([O2_pos, Ca1_pos, Ca2_pos, Ca3_pos, H2O_pos])
            x_pos = np.array([all_pos[k] for k in rangenames])

            # here is where we remembered the master line position for each feature so we can subtract it off and know how much we moved
            x_poses.append((x_pos.transpose() - x_pos[:, -1]).transpose())
            known_wave = x_pos[:,-1]


            for i in range(len(stars)-1):
                # create a new list wavelength solution that goes from the wavelength mosasurus gives to the new wavelength we have moved to
                fit = np.polyfit(x_pos[:,i], known_wave, 1)
                fit1d = np.poly1d(fit)
                # save the new wavelength array, as well as the coefficients so that we can re-create this fine-tuned wavelength solution later on in detrendersaurus
                newwavelength = stars[i]['wavelength']
                stars[i][width]['wavelength_adjusted'] = fit1d(newwavelength)
                stars[i][width]['stretch'] = fit[0]
                stars[i][width]['shift'] = fit[1]

                if makeplot:
                    color = self.starcolor(self.stars[i])
                    plt.cla()
                    plt.scatter(x_pos[:,i], known_wave - x_pos[:,i], color=color)
                    plt.plot(newwavelength, fit1d(newwavelength) - newwavelength, label=self.stars[i], color=color)
                    plt.xlabel('Original Wavelength (angstroms)')
                    plt.ylabel('New - Original (angstroms)')
                    dw = np.mean(known_wave - x_pos[:,i])
                    plt.ylim(- 10, 10)
                    plt.title(self.stars[i])

                    # make a directory
                    widthdir = os.path.join(self.starDirectories[i], 'wavelengthupdate_{}'.format(width))
                    mkdir(widthdir)
                    pltfilename = os.path.join(widthdir, 'wavelengthupdate_{}.pdf'.format(prefix))
                    plt.savefig(pltfilename)
                    self.speak('saved wavelength update plot to {}'.format(pltfilename))

        #a = input("how is {}?".format(width))

        for i in range(len(stars)-1):
            adir = self.starDirectories[i]
            newextractedpathname = os.path.join(adir, 'updatedwavelengthextracted_{}.npy'.format(prefix))
            #newextractedpathname = os.path.join(adir, 'extracted_{}.npy'.format(prefix))
            np.save(newextractedpathname, stars[i])
            self.speak('saved updated wavelengths in {}'.format(newextractedpathname))

  # ??????? (not sure if I'm using this) ### FIX ME ###
  def markBad(self):
      '''mark bad time-wavelength-star data points as bad'''
      satlimit = 150000
      faintlimit = 1

      for star in self.stars:
          self.speak('marking bad points for {}'.format(star))

          #if self.widths == None:
          #      self.widths = np.sort([np.float(x.split('_')[-1].split('px')[0])
          #                          for x in self.cubes['raw_counts'][star].keys()])

          for w in self.widths:
              self.widthkey = '{:04.1f}px'.format(w)
              self.speak('{} start ok'.format(np.sum(self.cubes['ok'][star][self.widthkey])))
              # mark wavelengths where the width is zero (or nearby)
              buffer = 10 # go this many pixels beyond borders
              undefined = self.cubes['width'][star][self.widthkey].sum(0) > 0
              bad = np.convolve(undefined, np.ones(buffer).astype(np.float)/buffer, mode='same')
              self.cubes['ok'][star][self.widthkey] *= (bad >= 0.9)

              # mark saturated values as not cool
              #self.cubes['ok'][star][self.widthkey] *= (self.cubes['peak'][star][self.widthkey] < satlimit)

              # mark things exceed the cosmic threshold as not cool

              #cosmics = self.squares['cosmicdiagnostic'][star][self.widthkey]
              #self.cubes['ok'][star][self.widthkey] *= cosmics[:,np.newaxis] < self.obs.cosmicAbandon
              # (DOUBLE CHECK THIS ISN'T A KLUDGE!)

              # mark things with really weird shifts as bad
              if self.shift:
                  shifts = self.squares['shift'][star][self.widthkey]
                  shifts -= np.median(shifts)
                  bad = np.abs(shifts) > 10*zachopy.oned.mad(shifts)
                  '''plt.figure('check shifts')
                  plt.cla()
                  plt.plot(shifts, color='gray')
                  plt.scatter(np.arange(len(shifts))[bad], shifts[bad], color='red')
                  plt.draw()'''
                  self.cubes['ok'][star][self.widthkey] *= (bad == False)[:,np.newaxis]
                  #self.input('plotting shifts for {}'.format(star))

              # mark
              self.cubes['centroid'][star][self.widthkey] -= np.median(self.cubes['centroid'][star][self.widthkey][self.cubes['ok'][star][self.widthkey]])
              self.speak('{} end ok'.format(np.sum(self.cubes['ok'][star][self.widthkey])))

  ### FIX ME ### (looks like this was made before multiple widths became a thing)
  def roughLC(self, target=None, comps=None, wavelengths=None, **kwargs):
      '''construct a rough LC, over a given wavelength bin'''

      # if comps is just one-element, make it into a list
      target = self.target
      comps = self.comparisons

      self.speak('updating the rough LC with target {0} and comps {1} and wavelength range {2}'.format(target, comps, wavelengths))

      w = self.spectral['wavelength']
      if wavelengths is None:
          wavelengths = [np.min(w), np.max(w)]
      blueenough = (w <= np.max(wavelengths))
      redenough = (w >= np.min(wavelengths))
      waveok = blueenough*redenough
      targetflux = self.cubes['raw_counts'][target][:,waveok].sum(-1)
      comparisonflux = self.cubes['raw_counts'][comps].sum(0)[:,waveok].sum(-1)

      self.speak('comp/targ is typically {0}'.format(np.median(comparisonflux/targetflux)))
      self.temporal['lc'] = targetflux/comparisonflux
      self.temporal['lc'] /= np.median(self.temporal['lc'])
      self.speak('rms is {0}'.format(np.std(self.temporal['lc'])))

      plt.ion()
      x, y = self.temporal['bjd'], self.temporal['lc']
      ok = self.temporal['ok']
      plt.plot(x[ok], y[ok], **kwargs)

  # ????????
  def doubleplot(self,
                    binsize=500,
                    wavemin=4500, wavemax=9000,
                    divide=False, ylim=None,
                    median=False,
                    title=None):

        self.populate()
        c = self
        name = c.obs.target.name
        date = c.obs.night.name
        wavelengths = np.arange(wavemin, wavemax, binsize)



        lcs = []

        import zachopy.cmaps
        blue, red = 'indigo', 'darkorange'
        cmap = zachopy.cmaps.one2another(blue, red)

        plt.ioff()
        for i, w in enumerate(wavelengths):
            c.roughLC(wavelengths=[w-binsize/2, w+binsize/2])
            ok = np.array(c.temporal['ok'])
            lc = c.temporal['bjd']+0.0, c.temporal['lc'] +0.0, ok
            lcs.append(lc)

        offset = int(c.temporal['bjd'][0])
        plt.figure(name, figsize=(8,6), dpi=70)
        if median:
            gs = plt.matplotlib.gridspec.GridSpec(3,1,
                        height_ratios=[1,1,.5], hspace=0.02)
        else:
            gs = plt.matplotlib.gridspec.GridSpec(2,1,
                        height_ratios=[1,.5], hspace=0.02)
        axlc = plt.subplot(gs[-1])

        plt.xlabel('BJD - {0}'.format(offset))
        plt.ylabel('Relative Flux')


        n, m = len(lc[2]), len(wavelengths)
        image, imageok = np.zeros((n,m)),  np.zeros((n,m))

        for i, lc in enumerate(lcs):
            fraction = i/(len(wavelengths) - 1.0)
            t = lc[0] - offset
            flux = lc[1]
            ok = lc[2]
            image[:,i] = flux + 0.0
            imageok[:,i] = ok + 0
            #print ok
            plt.plot(t[ok], flux[ok], alpha=binsize/500.0, linewidth=1,color=cmap(fraction))

        if ylim is None:
            valid = np.nonzero((imageok*np.isfinite(image)).flatten())[0]
            ylim = np.percentile(image.flatten()[valid], [1,99])

        plt.ylim(*ylim)
        plt.xlim(t.min(), t.max())


        if divide:
            image /= np.median(image, 1)[:,np.newaxis]

        axim = plt.subplot(gs[0])
        kw = dict(interpolation='nearest', cmap='gray',
                                vmin=ylim[0], vmax=ylim[1], aspect='auto',
                                extent=[min(t), max(t),
                                        (min(wavelengths) - binsize/2.0)/10, (max(wavelengths) + binsize/2.0)/10],
                                origin='lower')
        axim.imshow(image.T, **kw)


        plt.setp(axim.get_xticklabels(), visible=False)
        if title is None:
            title = '{name} with {instrument}\n[from {0}nm ({blue}) to {1}nm ({red}) in {2}nm-wide bins]'.format( wavemin/10, wavemax/10, binsize/10, name=name,blue=blue, red=red, instrument=c.obs.instrument.name)
        plt.title(title)
        plt.ylabel('Wavelength (nm)')

        if median:
            axmed = plt.subplot(gs[1])
            divided = (image/np.median(image, 1)[:,np.newaxis])
            kw['vmin'], kw['vmax'] = np.percentile(divided, [5,95])
            axmed.imshow(divided.T, **kw)
            plt.setp(axmed.get_xticklabels(), visible=False)

        plt.draw()

        if divide:
            filename = '{0}_binto{1}_{2}_normalized.pdf'.format(name,binsize, date)
        else:
            filename = '{0}_binto{1}_{2}.pdf'.format(name,binsize, date)
        plt.savefig(filename)

  # keep - displays the entire dataset, by making a movie of the individual spectra
  def movieCube(self,   fps=30, # how many frames per second
                        bitrate=1800*20, # bitrate (this seems to work well),
                        stride=10, **kw):# each frame will skip over this many timepoints):
      '''Create movie of the spectral cube.'''

      # set up a movie file, which will be populated with plots
      metadata = dict(artist='Z.K.B.-T.')
      self.writer = matplotlib.animation.FFMpegWriter(fps=fps, metadata=metadata, bitrate=bitrate)

      # set the filename based on whether cube was shifted
      if self.shift:
          filename = os.path.join(self.directory, 'shifted_cube_{0}stars_{1}spectra_{2}stride.mp4'.format(self.numberofstars, self.numberoftimes, stride))
      else:
          filename = os.path.join(self.directory, 'cube_{0}stars_{1}spectra_{2}stride.mp4'.format(self.numberofstars, self.numberoftimes, stride))


      if os.path.exists(filename):
          self.speak('a movie already exists at {}; skipping'.format(filename))
          return

      # plot the first spectrum, to set up the plots
      plt.ioff()
      try:
          # make sure we're restarting from scratch
          del self.ax_spectra
      except:
          pass
      self.plotSpectra(0, remake=True, **kw)

      # make the movie
      with self.writer.saving(self.figure, filename, self.figure.get_dpi()):
          # loop over exposures
          for i in np.arange(self.numberoftimes)[::stride]:
              self.speak('plotting {0}/{1}'.format(i, self.numberoftimes))
              self.plotSpectra(i)
              self.writer.grab_frame()

      # finish and display
      self.speak('saved movie to {0}'.format(filename))
      os.system('open {0}'.format(filename))

  # keep -- gives each star a unique color
  def starcolor(self, s):
    '''return a color for a particular star'''
    number = {s:i for i,s in enumerate(self.stellar['aperture'])}[s]
    return starcm(number/(self.numberofstars-1.0))


  # plots the spectra (and other things) for a single time-point (couples to movieCube)
  def plotSpectra(self, which, remake=False, wavelengthscale=10.0, figsize=None):
        '''For the ith component in the cube, plot all the spectra.'''

        # these (global) values will be plotted along the right side of the plot
        self.globallinekeys = ['airmass', 'rotatore']

        # these (star-by-star) values will be plotted along the right side
        self.starlinekeys = ['sky', 'width', 'centroid']#, 'shift']# 'cosmicdiagnostic', #, 'lc']
        if self.shift:
            self.starlinekeys = self.starlinekeys + ['shift']

        # these are the combination of linekeys
        self.linekeys = self.starlinekeys + self.globallinekeys

        # set up a plot figure
        #if figsize == None:
        #    figsize = (8 + self.numberofstars*4,24)
        figsize = (12, 8)
        self.figure = plt.figure('spectra', figsize=figsize, dpi=100)

        try:
            # determine whether or not the axes have been set up
            self.ax_spectra
            assert(remake == False)
        except (AssertionError, AttributeError):

            # create dictionaryies of axes and plots
            self.ax_spectra, self.ps = {}, {}

            # make gridspec structure with a row for each cubekey to plot, and a column for each star
            gs = plt.matplotlib.gridspec.GridSpec(len(self.cubekeys),self.numberofstars,hspace=0.08,wspace=0.2, left=0.1,right=0.75, bottom=0.1, top=0.9)

            # make a gridspec structure for vertical plots along the right side
            gstemporal = plt.matplotlib.gridspec.GridSpec(1, len(self.globallinekeys) + len(self.starlinekeys), left=0.78, right=0.98, wspace=0.05, bottom=0.1, top=0.9)
            sharex=None
            #ok = self.temporal['ok']

            # loop over the stars (columns)
            for istar, star in enumerate(self.stars):

                # figure out the self.widths

                #if self.widths == None:
                #    self.widths = np.sort([np.float(x.split('_')[-1].split('px')[0])
                #                          for x in self.cubes['raw_counts'][star].keys()])

                for i in range(len(self.cubekeys)):

                    # what's the key for the xaxis
                    k = self.cubekeys[i]


                    # add an axes for this key (row) and star (column)
                    try:
                        self.ax_spectra[k].append(plt.subplot(gs[i,istar], sharex=sharex))
                    except:
                        self.ax_spectra[k] = [plt.subplot(gs[i,istar], sharex=sharex)]
                    sharex = self.ax_spectra[k][0]
                    self.ax_spectra[k][0].set_ylabel(k)

                    # loop over self.widths, and plot them
                    w = self.width

                    # a mask setting which data points are okay
                    fine = self.cubes['ok'][star][which,:]

                    # pull out the y-values for this spectrum plot
                    thisy = self.cubes[k][star][self.cubes['ok'][star]]
                    thisy = thisy[np.isfinite(thisy)]
                    # make some adjustments, set ylim, and remove xlabels
                    if k in ['sky', 'peak', 'raw_counts']:
                        ylim = (0, np.percentile(thisy,99.9)*1.2)
                    elif k in ['centroid']:
                        ylim =  [-10.0, 10.0]#np.percentile(thisy,[1,99])
                    elif k in ['width']:
                        ylim = (0, 10.0)
                    self.ax_spectra[k][istar].set_ylim(*ylim )
                    plt.setp(self.ax_spectra[k][istar].get_xticklabels(), visible=False)

                # turn xlabels back on (just for bottom row)
                plt.setp(self.ax_spectra[k][istar].get_xticklabels(), visible=True)
                self.ax_spectra[k][istar].set_xlabel('Wavelength (nm)')

            # set the xlimits of the spectral plots
            self.ax_spectra[k][0].set_xlim(np.min(self.spectral['wavelength']/wavelengthscale), np.max(self.spectral['wavelength']/wavelengthscale))

            # now, plot the vertical plots along the right
            kw = dict(linewidth=1, alpha=1.0)
            sharey = None
            self.bars = []

            # do the global line keys
            for i, l in enumerate(self.linekeys):
                # create the appropriate axis
                ax = plt.subplot(gstemporal[i], sharey=sharey)
                sharey=ax
                if l in self.globallinekeys:
                    x = self.temporal[l]
                    ax.plot(x, np.arange(self.numberoftimes), color='black', **kw)
                    #ax.set_xlim(np.nanmin(self.temporal[l]),np.nanmax(self.temporal[l]))
                else:
                    for s in self.stars:
                        x = self.squares[l][s]
                        ok = np.isfinite(x)#self.cubes['ok'][s][:,:].max(1)
                        if l == 'cosmicdiagnostic':
                            x = np.log(x)
                        ax.plot(x[ok], np.arange(self.numberoftimes)[ok],
                                    color=self.starcolor(s),
                                    **kw)
                        print(x)
                        ax.set_xlim(np.nanmin(x[ok]), np.nanmax(x[ok]))

                # tidy up the plot
                ax.set_xlabel(l,rotation=45)
                self.bars.append(ax.axhline(which, color='gray', linewidth=4, alpha=0.5))
                plt.setp(ax.get_yticklabels(), visible=False)
                plt.setp(ax.get_xticklabels(), visible=False)
                ax.get_xaxis().get_major_formatter().set_useOffset(False)

                if i == int(len(self.linekeys)/2):
                    self.timestamp = ax.set_title('{0}: {1} {2}'.format(self.obs.target.name, self.temporal['ut-date'][which], self.temporal['ut-time'][which]))

            ax.set_ylim(self.numberoftimes, -1)
            for istar, star in enumerate(self.stars):
                self.ps[star] = {}
                for k in self.cubekeys:
                    try:
                        self.ps[star][k]
                    except KeyError:
                        self.ps[star][k] = self.ax_spectra[k][istar].plot(self.spectral['wavelength'], self.cubes[k][star][which,:], color=self.starcolor(star), **kw)[0]

        # set the position of the bars
        for b in self.bars:
            b.set_ydata(which)

        # set the spectra data for the plots
        for istar, s in enumerate(self.stars):
            for k in self.cubekeys:
                fine = self.cubes['ok'][s][which,:]
                x = self.spectral['wavelength'][fine]/wavelengthscale
                y = self.cubes[k][s][which,:][fine]
                self.ps[s][k].set_data(x, y)
                #self.speak('median {} value for star {} is {}'.format(k, s, np.median(y)))
            self.ax_spectra[self.cubekeys[0]][istar].set_title('image {0},\nstar {1}, aperture {2}'.format(self.temporal['fileprefix'][which], istar, s.replace('aperture_', '')))
        self.timestamp.set_text('{0}: {1} {2}'.format(self.obs.target.name, self.temporal['ut-date'][which], self.temporal['ut-time'][which]))

  ### FIX ME! ### (only partially switched back to single width)
  def shiftCube(self,plot=False):
    '''Shift all the spectra for a star to line up in wavelength.'''
    self.speak('shifting all spectra to line up at the calcium triplet')


    if plot:
        plt.figure('shifting spectra')
        plt.ion()

    # select a narrow wavelength range near the Ca triplet
    l, r = self.obs.instrument.extractiondefaults['correlationRange']
    l, r = 7500, 7750
    # define the range for this correlation
    left = np.argmin(np.abs(self.spectral['wavelength'] - l))
    right = np.argmin(np.abs(self.spectral['wavelength'] - r))
    width = self.obs.instrument.extractiondefaults['correlationSmooth']
    triplet = self.obs.instrument.extractiondefaults['correlationAnchors']

    c = self.cubes['raw_counts']

    # create a template to shift relative to (some fake, broadened absorption lines)
    wave = self.spectral['wavelength'][left:right]
    '''
    start = np.ones_like(wave).astype(np.float)
    for t in triplet:
      start *= 1 - np.exp(-(wave - t)**2/2/width**2)
    start -= np.median(start)
    '''

    def subtractcontinuum(y):
        return y - np.mean(y)

    # create a master template
    wave = self.spectral['wavelength'][left:right]
    masterexposure = 0
    master = self.cubes['raw_counts'][self.target][masterexposure, left:right]
    start = subtractcontinuum(master)

    # create an a structure to store the shifts in
    shifts = {}
    for star in self.stars:
        shifts[star] = np.zeros(self.numberoftimes)

    # calculate all the shifts
    for i in range(self.numberoftimes):
      for star in self.stars:

          # pull out the spectrum for one star,
          # at one extraction width, at one time point,
          # and trim it to a narrow range around the correlation anchors
          spectrum = c[star][i,left:right]

          # subtract its continuum
          this = subtractcontinuum(spectrum)

          # cross correlation with the anchors
          xc = np.correlate(this, start, 'same')

          # do a quadratic fit to estimate the peak
          x = np.arange(len(xc))
          coeff = np.polynomial.polynomial.polyfit(x[xc.argmax()-5:xc.argmax()+5], xc[xc.argmax()-5:xc.argmax()+5], 2)
          fit = np.polynomial.polynomial.Polynomial(coeff)
          der = fit.deriv()
          peak = der.roots()

          # calculate the offset from the peak
          offset = peak - len(wave)/2 + 1

          # now, use interpolation to shift all the arrays
          pixels = np.arange(self.numberofwavelengths)
          for key in self.cubes.keys():

              # should this be 'linear' or 'nearest'?
              interpolation = scipy.interpolate.interp1d(
                                pixels,self.cubes[key][star][i,:],
                                kind='linear',
                                bounds_error=False,
                                fill_value=0.0)

              interpolated = fluxconservingresample(pixels,
                                                    self.cubes[key][star][i,:],
                                                    pixels + offset,
                                                    treatnanas=0.0)

              self.cubes[key][star][i,:] = interpolated
              #interpolation(pixels + offset)
              #self.speak('shifting [{0}] by {1}A'.format(key, offset))

          if plot:
            # set up the plotting axes
            if i == 0:
              fi,ax = plt.subplots(2,1)
            else:
              for a in ax:
                a.cla()

            # plot the correlation function in the top panel
            ax[0].plot(x, xc, alpha=0.3)
            ax[0].scatter(peak, fit(peak))
            ax[0].set_xlim(0,len(x))

            # plot the spectrum
            ax[1].plot(wave, start/np.std(start), color='black')
            ax[1].plot(wave, start/np.std(start), color='black', alpha=0.2, linewidth=5)
            ax[1].plot(wave, this/np.std(this), color='gray', alpha=0.2, linewidth=5)
            new = subtractcontinuum(self.cubes['raw_counts'][star][i,left:right])
            ax[1].plot(wave, new/np.std(new), color='green', alpha=0.9, linewidth=2)
            ax[1].set_xlim(wave.min(), wave.max())
            ax[1].set_autoscaley_on
            ax[0].set_title('star {0}; {2}/{3} spectra'.format(star, len(self.cubes[key]), i+1, len(self.cubes[key][star][:,0])))
            ax[0].axvline(len(x)/2.0)
            #plt.draw()
            #self.input('?')
            pltdir = os.path.join(self.directory, star, 'shifts')
            mkdir(pltdir)
            pltfilename = os.path.join(pltdir, 'shift_{}.pdf'.format(self.obs.fileprefixes['science'][i]))
            plt.savefig(pltfilename)
          shifts[star][i] = offset
      self.speak( "shift = {4} for star {0}; {2}/{3} spectra".format(star, len(self.cubes[key]), i+1, len(self.cubes[key][star][:,0]), offset))

    self.squares['shift'] = shifts
    #self.temporal['shift'] = self.squares['shift'][self.obs.target[0],:]

  ### FIX ME (not used)
  def convolveCube(self,width=0.5):
    '''Take a spectral cube, convolve it with a Gaussian (NOT USED!).'''
    c = self.flux
    w = self.spectral['wavelength']
    plt.ion()
    x = np.arange(-width*10, width*10)
    gauss = lambda x: np.exp(-0.5*x**2/width**2)
    new_cube = np.copy(c)
    for star in range(len(c[:,0,0])):
        plt.figure('convolution')

        for i in range(len(c[0,:,0])):
            new_cube[star,i,:] = np.convolve(c[star,i,:], gauss(x),'same')
            plt.plot(c[star,i,:])
            plt.plot(new_cube[star,i,:])

    self.flux = new_cube

  def createBins(self, binsize=250, remake=False):
      '''Take a spectral cube, and creates big wavelength bins.'''
      self.populate()
      self.speak("Binning the spectral cube to {0}A binsize.".format(binsize))
      #self.shiftCube(plot=True)
      #self.show()
      wavelength = self.spectral['wavelength']
      #dnpixelsdw = self.spectral['dnpixelsdw']

      # the shape of the cube
      nStars, nTimes, nWaves = self.cubes['raw_counts'].shape
      plt.ion()

      self.binned_cubes = {}


      # define bins
      bin_starts = np.arange(wavelength.min(), wavelength.max() - binsize/2.0, binsize)
      bin_ends = bin_starts + binsize
      bin_centers = (bin_starts + bin_ends)/2.0
      nBins = bin_centers.size
      bin_ok = np.zeros((nStars, nBins))
      satlimit = 150000
      faintlimit = 2500

      # loop through stars, determining which bins are okay for each
      for k in self.cubes.keys():
          self.speak('binning {0} to {1}A bins'.format(k, binsize))
          shape = (self.numberofstars, self.numberoftimes, self.numberofwavelengths/binsize, binsize)

          # integrate spectrum over bin (dndw factor accounts for the fractional pixels represented by resampled cube spectrum
          # !!!!! DO I NEED TO INCLUDE DNPIXELSDW??
          # this was the old version:
          #                      self.binned_cubes[newkey][star, time, bin] = scipy.integrate.trapz(wholespectrum[mask]*dnpixelsdw[mask], wavelength[mask])
          #if k == 'ok':
          #      dndw = 1
          #else:
          #  dndw = self.spectral['dnpixelsdw'].reshape(1, 1, self.numberofwavelengths)
          self.binned_cubes[k] = (self.cubes[k]*dndw).reshape(shape).sum(-1)
          if k=='width' or k =='centroid' or k =='peak':
              self.binned_cubes[k] /= (dndw).reshape((1,1, self.numberofwavelengths/binsize, binsize)).sum(-1)
          if k=='ok':
              self.binned_cubes[k] = (self.binned_cubes[k]/binsize).astype(np.bool)
      #
      self.bin_centers = bin_centers
      self.binned_cubes = astropy.table.Table(self.binned_cubes)
      #np.save(binned_filename, (self.binned_cubes, self.bin_centers, self.binned_cubes['ok']))
      #self.speak("    Saved binned cube to {0}".format(binned_filename))

  def correctBins(self, **kw):
    '''Use comparison stars to correct for atmospheric losses, create a self.binned_corrected.'''

    self.populate()
    self.speak('using comparison stars {0} to correct for atmospheric losses'.format(self.comparisons))
    wavelength = self.bin_centers
    vmin = 0.98
    vmax = 1.02
    nStars, nTimes, nWaves = self.binned_cubes['raw_counts'].shape

    # create empty correction and uncertainty arrays
    correction, uncertainty = np.ones((nTimes, nWaves)), np.ones((nTimes, nWaves))

    def weightedsum(array):
        return ((array*self.binned_cubes['ok'])[self.comparisons,:,:].sum(0)/(self.binned_cubes['ok'])[self.comparisons,:,:].sum(0))

    correction = weightedsum(self.binned_cubes['raw_counts'])
    uncertainty = np.sqrt(weightedsum(self.binned_cubes['raw_counts'] + self.binned_cubes['sky']))

    self.binned_correction = np.ma.MaskedArray(correction, mask=((self.binned_cubes['ok']==False).sum(0).astype(np.bool)), fill_value=np.nan)
    self.binned_correction_uncertainty = np.ma.MaskedArray(uncertainty, mask=((self.binned_cubes['ok']==False).sum(0).astype(np.bool)), fill_value=np.nan)

    # normalize the correction spectrum to be close to one
    mediancompositespectrum = np.ma.median(self.binned_correction, 0)
    self.binned_correction /= mediancompositespectrum.reshape(1,nWaves)
    self.binned_correction_uncertainty /= mediancompositespectrum.reshape(1,nWaves)
    #self.display.one(self.binned_correction.filled(), clobber=True)
    #self.display.one(self.binned_correction_uncertainty.filled())
    self.binned_cubes['corrected'] = self.binned_cubes['raw_counts']/self.binned_correction
    photonnoise = np.sqrt(self.binned_cubes['raw_counts'] + self.binned_cubes['sky'])/self.binned_cubes['raw_counts']
    correctionnoise = self.binned_correction_uncertainty
    self.binned_cubes['uncertainty'] = np.sqrt(photonnoise**2 + correctionnoise**2)

  def imageBins(self, **kw):
      self.createBins(**kw)
      self.correctBins(**kw)

      figure = plt.figure(figsize=(10,10), dpi=70)
      gs = plt.matplotlib.gridspec.GridSpec(self.numberofstars, 2)

      kw = dict(cmap='gray')
      ax=None
      for i in range(self.numberofstars):
          ax = plt.subplot(gs[i,0], sharex=ax, sharey=ax)
          ax.imshow(self.binned_cubes['raw_counts'][i], **kw)
          ax = plt.subplot(gs[i,1], sharex=ax, sharey=ax)
          ax.imshow(self.binned_cubes['corrected'][i], **kw)

  def makeMeanSpectrum(self, plot=False):
    self.populate()
    wavelength = self.spectral['wavelength']
    spectrum = np.median(self.cubes['raw_counts'][self.obs.target,:,:],1).flatten()
    assert(len(spectrum) == len(wavelength))
    if plot:
        fi, ax = plt.subplots(1)
        unit = 10.0
        ax.plot(wavelength/unit, spectrum*unit)
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Flux (photons/nm/exposure)')
    self.speak("saving median spectrum to")
    filename = os.path.join(self.directory, 'medianSpectrum.npy')
    self.speak(filename)
    np.save(filename, (wavelength, spectrum))

  def makeLCs(self,binsize=250, remake=False):
    '''Wrapper to go from extracted spectra to binned, multiwavelength lightcurves.'''

    self.populate()
    # make (and save) and mean spectrum, for plotting comparisons at later steps
    self.makeMeanSpectrum()

    # pick the target star
    target = self.target
    comparisons = self.obs.comparisons

    # bin the cube into manageable wavelength bins
    self.createBins(binsize=binsize)

    # use comparison star(s) to divide out flux losses
    self.correctBins()

    # setup
    nStars, nTimes, nWaves = self.binned_cubes['corrected'].shape
    bw = self.bin_centers
    binsize = bw[1] - bw[0]
    bin_starts = bw - binsize/2.0
    bin_ends = bw + binsize/2.0


    lcDirectory = os.path.join(self.directory,  "chromatic{binsize:05.0f}/".format(binsize=binsize))
    mkdir(lcDirectory)
    lcDirectory = os.path.join(lcDirectory, 'originalLCs/')
    mkdir(lcDirectory)

    self.lcs = []

    # loop through wavelength bins
    for wave in range(nWaves):
      left, right = bin_starts[wave], bin_ends[wave]
      lcfilename =  os.path.join(lcDirectory, '/{0:05d}to{1:05d}.lightcurve'.format(np.int(left), np.int(right)))

      # is there *any* good data at this wavelength?
      if self.binned_cubes['ok'][target,:,wave].any():

        # determine where the target star is ok and the correction is ok (on a time-by-time basis)
        ok = ((self.binned_cubes['ok'][target,:,wave] != False).flatten()*(self.binned_correction.mask[:,wave] == False).flatten())

        # report what's being made
        self.speak('making light curve for {0} to {1}A, with {2} good points'.format(bin_starts[wave], bin_ends[wave], np.sum(ok != False)))

        # create an empty LC object
        lc = astropy.table.Table()
        lc['bjd'] = self.temporal['bjd']
        lc['flux'] = self.binned_cubes['corrected'][target,:,wave].flatten()/np.median(self.binned_cubes['raw_counts'][target,:,wave].flatten())
        lc['uncertainty'] = self.binned_cubes['uncertainty'][target,:,wave].flatten()
        lc['ok'] = ok.astype(np.int)

        # pull out global values
        for key in ['airmass', 'rotatore']:
            lc['{0}'.format(key)] = self.temporal[key]

        # pull out the star-by-star (wavelength-independent) quantities
        for key in ['width', 'centroid', 'shift']:
            try:
                lc['{0}_target'.format(key)] = self.squares[key][target]
                for comparison in comparisons:
                    lc['{0}_star{1:02.0f}'.format(key, comparison)] = self.squares[key][comparison]
            except KeyError:
                self.speak("{} couldn't be found!".format(key))
        # pull out the star-by-star wavelength specific values
        for key in ['sky', 'peak']:
            lc['{0}_target'.format(key)] = self.binned_cubes[key][target,:,wave]
            for comparison in comparisons:
                lc['{0}_star{1:02.0f}'.format(key, comparison)] = self.binned_cubes[key][comparison,:,wave]

        # pull out the star-by-star wavelength specific values that should be measured relative to the more global values
        for key in ['width', 'centroid']:
            lc['d{0}_target'.format(key)] = self.binned_cubes[key][target,:,wave] - lc['{0}_target'.format(key)]
            for comparison in comparisons:
                lc['d{0}_star{1:02.0f}'.format(key, comparison)] = self.binned_cubes[key][comparison,:,wave] - lc['{0}_star{1:02.0f}'.format(key, comparison)]

        #lc.populate(bjd, flux, uncertainty, **lc)
        table = astropy.table.Table(lc)
        table['bjd'].format = '.10f'

        #table = table[table['ok'].astype(np.bool)]
        # REMOVED TO MAKE SURE MASKING IS EASIER AT LATER STEP


        table.write(lcfilename, format='ascii.fixed_width', bookend=False)
        self.speak('saved light curve to')
        self.speak('{0}'.format(lcfilename))

        '''for key in self.binned_cubes.keys():
        if key != 'corrected' and key != 'error':
          dict[key+'_target'] = self.binned_cubes[key][target,ok,wave].flatten()
          for comparison in self.obs.comparisons:
            dict[key+'_comparison{0:02.0f}'.format(comparison)] = self.binned_cubes[key][comparison,ok,wave].flatten()
        for k in keystoinclude:
          if k == 'ok':
              continue
          if k == 'airmass':
              newkey = k
          else:
              newkey = 'global_{0}'.format(k)
          dict[k] = self.temporal[k][ok]

          assert(np.isfinite(flux[ok]).all())
          assert(np.sum(ok) > 0)
          self.speak(bjd.flatten()[ok].size)
          lc.populate(bjd[ok], flux[ok], error[ok], **dict)

          #lc.plot()
          self.speak(lc)
          lc.save()
          self.lcs.append(lc)'''

  def loadLCs(self, binsize=250):
    lcDirectory = os.path.join(self.directory, 'lc_binby' + ('%d' % binsize) + '/')
    g = glob.glob(os.path.join(lcDirectory, 'lc_*.npy'))
    wavelengths = []
    lcs = []
    for file in g:
      lc = LC(self.obs, filename=file)
      if lc.lc is not None:
        wavelengths.append(lc.wavelength)
        lcs.append(lc)
    self.lcs = np.array(lcs)
    self.lcs = self.lcs[np.argsort(wavelengths)]
    return self.lcs




  def imageTarget(self, title=None, vmin=None, vmax=None):
    '''Make an image of the input cube.'''
    if title is None:
      title = self.obs.target.name + ' | ' + self.obs.night
    self.speak("   Trying to image bins for " + title)
    bin_centers = self.bin_centers
    bin_ok = self.binned_cubes['ok']
    target = self.target_raw
    '''Take binned lightcurves, image them. Returns the bins with median correction subtracted.'''
    plt.ion()
    if vmin == None:
      vmin = 0.985
    if vmax == None:
      vmax = 1.005
    nTimes, nWaves = target.shape
    fi, ax = plt.subplots(4,1, figsize=(4.5,12), sharex=True, sharey=True)
    plt.subplots_adjust()
    ax[0].set_title(title)
    binsize = bin_centers[1] - bin_centers[0]
    self.bin_starts = bin_centers - binsize/2.0
    self.bin_ends = bin_centers + binsize/2.0
    nBins = bin_centers.size

    normalized = np.ones_like(self.target_raw)
    targets = [self.target_raw, self.target_corrected, self.target_median, self.target_divided]
    names = ['Raw Detected Flux', 'Divided by Comparison', 'Median Transit', 'Divided by Transit']
    for i in range(len(targets)):
      target = targets[i]

      medianspectrum = np.median(target,0)*np.ones_like(target)
      normalized = target/medianspectrum
      ax[i].imshow(normalized, vmin=vmin, vmax=vmax, extent=[self.bin_starts.min(), self.bin_ends.max(), 0, nTimes], aspect='auto', cmap='gray', interpolation='nearest')
      ax[i].set_ylabel(names[i])
    ax[-1].set_xlabel('Wavelength (angstroms)')
    ax[-1].set_xlim(self.bin_starts.min(), self.bin_ends.max())
    for bin in range(nBins):
      if bin_ok[self.obs.target,bin] == 0:
        for a in ax:
          a.axvspan(self.bin_starts[bin], self.bin_ends[bin], alpha=0.7, color='white', edgecolor=None)
    plt.tight_layout(h_pad=0.0)
    plt.savefig(os.path.join(self.obs.reducer.extractionDirectory, + 'divided_{0}.pdf'.format(self.obs.night)))
    a = raw_input('imaging  target')

  def help(self):
      self.speak(cubehelp)


class LC():
  def __init__(self, obs, left=None, right=None, filename=None):
    self.obs = obs
    if filename is not None:
      self.load(filename)
    else:
      self.left = left
      self.right = right
      self.setup()


  def setup(self):
    self.wavelength = (self.left + self.right)/2.0
    self.binsize = self.right - self.left
    self.filename = os.path.join(self.directory, 'lc_binby' + ('%d' % self.binsize) + '/lc_{0:05d}to{1:05d}.npy'.format(np.int(self.left), np.int(self.right)))

  def populate(self, bjd, flux, error, **kwargs):
    # set up the column names for the light curve record array
    types = [('bjd', np.float), ('flux', np.float), ('error', np.float)]
    for key in kwargs.keys():
      types.append((key,np.float))

    # populate the columns with data
    self.lc = np.zeros(bjd.size, types)
    self.lc['bjd'] = bjd
    self.lc['flux'] = flux
    self.lc['error'] = error
    for key in kwargs.keys():
      self.lc[key] = kwargs[key]



  def save(self):
    np.save(self.filename, self.lc)

  def load(self, filename):
    self.left = np.float(filename.split('lc_')[-1].split('to')[-2])
    self.right = np.float(filename.split('.npy')[-2].split('to')[-1])
    self.setup()
    assert(self.filename == filename)
    self.lc = np.load(self.filename)


  def plot(self):
    try:
      for a in self.ax.values():
        a.cla()
    except:
      self.cubekeystoplot = ['flux', 'raw_counts_target', 'sky_target', 'airmass', 'width_target', 'centroid_target', 'peak_target']
      gs = plt.matplotlib.gridspec.GridSpec(len(self.cubekeystoplot), 1, height_ratios=[5,1,1,1,1,1,1], wspace=0, hspace=0)
      self.ax = {}
      for i in range(len(self.cubekeystoplot)):
        key = self.cubekeystoplot[i]
        try:
          sharex = self.ax[self.cubekeystoplot[0]]
        except:
          sharex = None
        self.ax[key] = plt.subplot(gs[i], sharex=sharex)
    for key in self.cubekeystoplot:
      self.ax[key].plot(self.lc['bjd'], self.lc[key], markersize=1, marker='o', alpha=0.3, color='black', linewidth=0)
      self.ax[key].set_xlim(self.lc['bjd'].min(), self.lc['bjd'].max())
      self.ax[key].set_ylabel(key)

    plt.draw()
    bla = raw_input(self.filename + '?')
