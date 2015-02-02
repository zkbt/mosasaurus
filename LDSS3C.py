# required things to add:
# estimate gain for amplifiers, convert everything to electrons early on


from Observation import Observation
from Night import Night
from Tools import *
from Display import Display
from CCD import CCD
from Calibration import Calibration
from Mask import Mask
from Aperture import Aperture



def reduce(filename= 'wasp94_140801.obs'):
	# generate an observation object
	obs = Observation(filename)
	# load the headers for this observation
	obs.loadHeaders()
	# create a night object, associated with this observation
	night = Night(obs)
	# create an observing log from all the image headers from this night
	night.obsLog()

	# set up the calibration
	calib = Calibration(obs)

	# loop through the CCD's needed for this observation, and make sure they are stitched
	#for n in obs.nNeeded:
	#		ccd = CCD(obs,n=n,calib=calib)
	#	ccd.createStitched(visualize=True)


	mask = Mask(calib)
	for a in mask.apertures:
		a.displayStamps(a.images)
		a.extractAll(remake=True)

def show(filename = 'wasp94_140801.obs', binsize=25, vmin=0.98, vmax=1.02, remake=False):
	obs = Observation(filename, nods9=True)
	cube = Cube(obs, remake=remake)
	cube.makeMeanSpectrum()

	#cube.loadSpectra()
	#cube.convolveCube(width=5.0)
	#cube.shiftCube()
	cube.makeLCs(binsize=binsize)

def divide(filename = 'wasp94_140801.obs', binsize=25, remake=False):
	obs = Observation(filename, nods9=True)
	cube = Cube(obs)
	cube.createBins(binsize=binsize, remake=remake)
	cube.correctBins()
	cube.divide()

def fit(filename = 'wasp94_140801.obs', binsize=25, remake=False):
	obs = Observation(filename, nods9=True)
	cube = Cube(obs)
	if remake:
		cube.makeLCs(binsize=binsize)
	cube.loadLCs(binsize=binsize)
	tm = fitTLC.TM(e=0.0, a_over_rs=0.055*1.496e13/1.45/6.96e10, b=0.17, u1=0.1, u2 = 0.1, rp_over_rs=np.sqrt(0.01197), w=np.pi/2.0, period=3.9501907, t0=2456416.40138)
	for lc in cube.lcs:
		tlc = fitTLC.TLC(lc.lc['bjd'], lc.lc['raw_counts'])
		tm.plot(tlc)
	return cube
