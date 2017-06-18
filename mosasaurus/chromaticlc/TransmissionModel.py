import glob, os, string, copy
from  zachopy.borrowed.mpfit.mpfit import mpfit


import astropy.io.fits, astropy.io.ascii, astropy.time
#from astropy.io import fits, ascii
#from astropy.time import Time

import matplotlib.pyplot as plt, numpy as np, matplotlib.animation
import zachopy.units as u
from zachopy.Talker import Talker
import zachopy.oned
tablefile = '/Users/zkbt/Dropbox/poptrans/data/exoplanets.org.txt'




class TransmissionModel(Talker):
    def __init__(self, planet='WASP94Ab'):
        Talker.__init__(self)
        data = astropy.io.ascii.read(tablefile)
        match = (data['NAME'] == planet).nonzero()[0]
        assert(len(match == 1))
        self.row = data[match]
        self.mass = self.row['MASS'][0] # in Jupiter
        self.radius = self.row['R'][0] # in Jupiter
        self.period = self.row['PER'][0]
        self.a_over_r = self.row['AR'][0]
        self.stellar_teff = self.row['TEFF']
        self.stellar_mass = self.row['MSTAR']
        self.stellar_radius = self.row['RSTAR']

        self.teq =  self.stellar_teff/np.sqrt(2*self.a_over_r)
        self.mu = 2.2
        self.surfacegravity = u.G*self.mass*u.Mearth/(self.radius*u.Rearth)**2
        self.scaleheight = u.k_B*self.teq/u.mp/self.mu/self.surfacegravity

        print self.__dict__


    def propagate(self):
        optical = (self.wavelength < 7000)*(self.wavelength > 4000)
        self.original_deltaradiusinrj = self.original_radiusinrj - np.min(self.original_radiusinrj[optical])
        self.scale_factor = self.scaleheight/self.original_scaleheight

        self.nscaleheights = self.original_deltaradiusinrj*u.Rjupiter/self.original_scaleheight
        self.deltaradiusinrj = self.original_deltaradiusinrj*self.scale_factor
        self.deltadepth = 2*self.radius*u.Rearth*self.deltaradiusinrj*u.Rjupiter/(self.stellar_radius*u.Rsun)**2
        self.depth = (self.radius*u.Rearth/self.stellar_radius/u.Rsun)**2

    def binas(self, spectrum):
        smoothed, wavelengths = [], []
        smoothed_nscaleheights = []
        for i in range(len(spectrum.wavelength)):
            # KLUDGE! YOU BETTER MAKE SURE YOUR UNITS MATCH UP!
            ok = (self.wavelength/10.0 >= spectrum.left[i])*(self.wavelength/10.0 < spectrum.right[i])
            smoothed.append(np.mean(self.deltadepth[ok]) + self.depth)
            smoothed_nscaleheights.append(np.mean(self.nscaleheights[ok]))

        self.smoothed_depth = np.array(smoothed)
        self.smoothed_wavelength = np.array(spectrum.wavelength)
        self.smoothed_nscaleheights = np.array(smoothed_nscaleheights)
        #assert(np.isfinite(self.smoothed_nscaleheights).all())
        return self.smoothed_nscaleheights

    def binned(self, binwidth=100):
        return zachopy.oned.binto(self.wavelength, self.nscaleheights, binwidth)

    def deviates(self, p, fjac=None, spectrum=None):
        self.nudge = p[0]
        model = (self.binas(spectrum) + self.nudge)*self.scaleheight

        dev = (spectrum.dradius - model)/spectrum.dradius_error
        dev = dev[np.isfinite(model)]
        return [0, dev]

    def chisq(self, spectrum):
        p0 = [0.0]
        parinfo = [{'value':0., 'fixed':0, 'limited':[1,1], 'limits':[-10,10], 'step':0.1}]
        self.fitted = mpfit(self.deviates, p0, parinfo=parinfo, functkw=dict(spectrum=spectrum))
        dev = self.deviates(self.fitted.params, spectrum=spectrum)[-1]
        chisq = np.sum(dev**2)
        dof = self.fitted.dof
        return chisq, dof



class Fortney(TransmissionModel):
    '''Can't be run on its own, needs filename defined in inherited classes.'''
    def __init__(self, planet='WASP94Ab', modelpath = '/Users/zkbt/Dropbox/fortneyhotjupiter/'):

        TransmissionModel.__init__(self, planet=planet)

        self.original_scaleheight = u.k_B*self.original_teq/u.mp/self.original_mu/self.original_surfacegravity
        self.data = astropy.io.ascii.read(self.filename)

        sort = np.argsort(self.data['wavelength'])
        self.data = self.data[sort]
        self.data = self.data[self.data['wavelength'] < 1.0]
        self.wavelength = self.data['wavelength']*10000.0 # to get it in angstrom
        self.original_radiusinrj = self.data['radiusinkm']*1e5/u.Rjupiter
        self.propagate()




class WithoutTiO(Fortney):
    def __init__(self, planet='WASP94Ab', modelpath = '/Users/zkbt/Dropbox/fortneyhotjupiter/'):
        # set up which file to use
        self.filename = modelpath + 'lambda_1500K_g10_noTiOVO.dat'
        self.color = 'DarkOrange'
        self.original_teq = 1500.0
        self.original_surfacegravity = 10.0*100.0
        self.original_mu = 2.2
        self.name = 'cloud-free model'
        # initialize the Fortney model
        Fortney.__init__(self, planet=planet)


class WithTiO(Fortney):
    def __init__(self, planet='WASP94Ab', modelpath = '/Users/zkbt/Dropbox/fortneyhotjupiter/'):
        # set up which file to use
        self.filename = modelpath + 'lambda_1500K_g10_wTiOVO.dat'
        self.color='RoyalBlue'
        self.original_teq = 1500.0
        self.original_surfacegravity = 10.0*100.0
        self.original_mu = 2.2
        self.name = 'cloud-free + TiO model'

        # initialize the Fortney model
        Fortney.__init__(self, planet=planet)


class Morley(TransmissionModel):
    def __init__(self, planet='WASP94Ab', modelpath = '/Users/zkbt/Dropbox/morleywasp94/'):


        TransmissionModel.__init__(self, planet=planet)

        self.original_scaleheight = 154138541.634 #self.scaleheight#u.k*self.original_teq/u.mp/self.original_mu/self.original_surfacegravity
        self.data = astropy.io.ascii.read(self.filename)

        sort = np.argsort(self.data['col1'])
        self.data = self.data[sort]
        self.data = self.data[self.data['col1'] < 1.0]
        self.wavelength = self.data['col1']*10000.0 # to get it in angstrom
        self.original_radiusinrj = self.data['col2']*1e5/u.Rjupiter
        self.propagate()

class NoClouds(Morley):
    def __init__(self, planet='WASP94Ab', modelpath = '/Users/zkbt/Dropbox/morleywasp94/'):
        # set up which file to use
        self.filename = modelpath + 'wasp94Ab-m0.0-nc.out'
        self.color = 'DarkOrange'
        self.name = 'cloud-free'

        # initialize the Fortney model
        Morley.__init__(self, planet=planet)

class fsed05(Morley):
    def __init__(self, planet='WASP94Ab', modelpath = '/Users/zkbt/Dropbox/morleywasp94/'):
        # set up which file to use
        self.filename = modelpath + 'wasp94Ab-m0.0-f0.5_100.out'
        self.color = 'DarkGreen'
        self.name = 'fsed=0.5'

        # initialize the Fortney model
        Morley.__init__(self, planet=planet)


class fsed20(Morley):
    def __init__(self, planet='WASP94Ab', modelpath = '/Users/zkbt/Dropbox/morleywasp94/'):
        # set up which file to use
        self.filename = modelpath + 'wasp94Ab-m0.0-f2_100.out'
        self.color = 'SeaGreen'
        self.name = 'fsed=2'

        # initialize the Fortney model
        Morley.__init__(self, planet=planet)

class fsed30(Morley):
    def __init__(self, planet='WASP94Ab', modelpath = '/Users/zkbt/Dropbox/morleywasp94/'):
        # set up which file to use
        self.filename = modelpath + 'wasp94Ab-m0.0-f3_100.out'
        self.color = 'MediumSeaGreen'
        self.name = 'fsed=3'
        # initialize the Fortney model
        Morley.__init__(self, planet=planet)

class Flat(Morley):
    def __init__(self, planet='WASP94Ab', modelpath = '/Users/zkbt/Dropbox/morleywasp94/'):
        # set up which file to use
        self.filename = modelpath + 'wasp94Ab-m0.0-f3_100.out'
        self.color = 'Gray'
        self.name = 'flat model'
        # initialize the Fortney model
        Morley.__init__(self, planet=planet)
        self.original_radiusinrj[:] = np.mean(self.original_radiusinrj[:])
        self.propagate()
