#!/usr/bin/env python
# relies on Jonathan Irwin's library for astronomy
# (which is amazing)
from astropy.coordinates import EarthLocation
import sys
import math
import lfa

def jdutc2bjdtdb(jd_utc, ra, dec, observatory='lco',
            pmra=0.0, pmde=0.0, plx=0.0, catep=2000.0, vrad=0):
    '''wrapper to use Jonathan Irwin's lfa code to calculate BJD_TDB
        inputs:
            jd_utc = time in JD_UTC
            ra = Right Ascension in deg
            dec = Dec in degrees
            observatory = where on Earth (as input to astropy EarthLocation)
            (optional)
            pmra, pmde = proper motions in projected arcsec/yr
            plx = parallax in arcsec
            catep = the epoch of the catalog position
            vrad = radial velocity

        neglects many important factors
        that would be important for getting
        below 1 second precision -- beware!'''

    # pull out the observatory parameters
    print "observatory is {0}".format('lco')
    #o = EarthLocation.of_site('lco')
    longitude = 0# o.longitude.rad
    latitude = 0#o.latitude.rad
    height = 0#o.height.to('m').value
    # (taken J's values for CTIO)
    temperat   = 283.16
    humidity   = 0.5
    pressure   = 777.1
    wavelength = 0.7

    utc = jd_utc - 2400000.5

    catra = ra * lfa.DEG_TO_RAD
    catde = dec * lfa.DEG_TO_RAD

    # Jonathan's observer structure
    obs = lfa.observer(longitude, latitude, height)

    # Atmospheric refraction.
    obs.refract(temperat, humidity, pressure, wavelength)

    # Figure out TT from given UTC.
    iutc = int(utc)
    ttmutc = obs.dtai(iutc, utc-iutc) + lfa.DTT;

    # Compute time-dependent quantities.
    obs.update(utc, ttmutc, lfa.OBSERVER_UPDATE_ALL);

    # Figure out total clock correction TDB-UTC.
    dclock = ttmutc + obs.dtdb;

    # New source structure.
    src = lfa.source(catra, catde, pmra, pmde, plx, vrad, catep);

    # Compute current BCRS position.
    (s, dsdt, pr) = obs.place(src, lfa.TR_MOTION);

    # Barycentric delay.
    delay = obs.bary_delay(s, pr);

    print "DELAY = Bary {0:.10f}s Clock {1:.10f}s Total {2:.10f}s".format(delay, dclock, delay+dclock)

    print "modified BJD(TDB) = {0:.10f}".format(utc+(delay+dclock)/lfa.DAY)
    return jd_utc + (delay+dclock)/lfa.DAY
