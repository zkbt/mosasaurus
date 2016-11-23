'''Tools for converting from times to BJD.'''

from astropy import time, coordinates as coord, units as u

def toBJD(  times, # astropy time object, with location + scale specified
            ra, dec, # ra and dec, in degrees
            verbose=True):
    '''get BJD from astropy times (and coordinates and locations)'''

    sky = coord.SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')
    ltt_bary = times.light_travel_time(sky)

    times_bary = times.tdb + ltt_bary
    if verbose:
        print """
            For JD_UTC {}
            at position {},
            from site {},
            BJD_TDB - JD_UTC is {} (days),
            and BJD_TDB is {}.
        """.format(   times.utc.jd,
                            sky.to_string('hmsdms'),
                            times.location.geodetic,
                            (times_bary.tdb.jd - times.utc.jd),
                            times_bary.tdb.jd)

    return times_bary

def test():
    '''make sure the toBJD function works.'''

    loc = coord.EarthLocation.from_geodetic(0, 0, 0)
    t = time.Time(2457000.0, format='jd', scale='utc', location=loc)
    bjd = toBJD(t, 0, 0)

    print """
    It might be a good idea to test this against
    http://astroutils.astronomy.ohio-state.edu/time/utc2bjd.html
    """
