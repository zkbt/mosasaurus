'''Tools for converting from times to BJD.'''

from astropy import time, coordinates as coord, units as u

def toBJD(
            times, # astropy time object, with location + scale specified
            skyposition, # astropy SkyCoord object, with RA + Dec in icrs
            verbose=True, # should we talk about the result?
         ):
    '''get BJD from astropy times (and coordinates and locations)'''

    # determine the light travel time to the barycenter
    ltt_bary = times.light_travel_time(skyposition)

    # convert into the TDB scale, and add the light travel time
    times_bary = times.tdb + ltt_bary

    if verbose:
        # print out the results, if you're curious
        print("""
        For JD_UTC {}
        at position {},
        from site {},
        BJD_TDB - JD_UTC is {:.6} (days) = {:.6} (minutes),
        and BJD_TDB is {}.
        """.format( times.utc.jd,
                    skyposition.to_string('hmsdms'),
                    times.location.geodetic,
                    (times_bary.tdb.jd - times.utc.jd),
                    (times_bary.tdb.jd - times.utc.jd)*24*60,
                    times_bary.tdb.jd))

    # return the times, as astropy times
    return times_bary


def visualizeLTT():
    '''
    make a cartoon drawing the Earth, the Sun, the direction to the star
    and indicating the light travel time to the barycenter
    '''

    # (not yet implemented)
    return "trombone noise"

def test():
    '''make sure the toBJD function works.'''

    # specify a location relative to Earth
    loc = coord.EarthLocation.from_geodetic(0, 0, 0)
    # specify a JD_UTC
    t = time.Time(2457000.0, format='jd', scale='utc', location=loc)
    # specify a location on the sky
    position = coord.SkyCoord(0.0*u.deg, 0.0*u.deg, frame='icrs')
    # calculate the bjd for this time
    bjd = toBJD(t, position, verbose=True)

    print("""
    I very strongly suggest testing this result against Jason Eastman's tool at
        http://astroutils.astronomy.ohio-state.edu/time/utc2bjd.html
    """)

if __name__ == '__main__':
    test()
