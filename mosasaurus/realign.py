
 for i in range(subny): # subny == number of pixels in disp. direction
        trdata[i]  = subdata [i,itrace[i]-spec_width:itrace[i]+spec_width]
        trmask[i]  = submask [i,itrace[i]-spec_width:itrace[i]+spec_width]
        trvar [i]  = variance[i,itrace[i]-spec_width:itrace[i]+spec_width]
        trbg  [i]  = bg      [i,itrace[i]-spec_width:itrace[i]+spec_width]

    trdata  = np.zeros((subny,spec_width*2))
    trmask  = np.zeros((subny,spec_width*2),dtype=bool)
    trvar   = np.zeros((subny,spec_width*2))
    trbg    = np.zeros((subny,spec_width*2))
    spldata = spi.RectBivariateSpline(range(subny),range(subnx),subdata, kx=3, ky=3, s=0)
    splmask = spi.RectBivariateSpline(range(subny),range(subnx),submask, kx=3, ky=3, s=0)
    splvar  = spi.RectBivariateSpline(range(subny),range(subnx),variance, kx=3, ky=3, s=0)
    splbg   = spi.RectBivariateSpline(range(subny),range(subnx),bg, kx=3, ky=3, s=0)
    for i in range(subny):
        xi          = (np.arange(trace[i]-spec_width,trace[i]+spec_width)+0.5)[:spec_width*2]
        yi          = np.zeros(len(xi))+i
        trdata[i]   = spldata.ev(yi,xi)
        trmask[i]   = np.round(splmask.ev(yi,xi),0).astype(bool)
        trvar [i]   = splvar .ev(yi,xi)
        trbg  [i]   = splbg  .ev(yi,xi) 


# variables I do not know:
#   I assume subdata, submask, and bg are the same as the parameters of the same names that I feed into the optimized spectrum
#   what are variance, itrace, trace, and spec_width?
#   you have a comment saying subny == number of pixels in dispersion direction, so is subnx the number of pixels in the cross-dispersion direction? if so, how is this different from spec_width?
#   it seems that spi.RectBivariateSpline is a standard scipy package, so I'm assuming that won't be a problem
#   and finally, I think I can figure this out if I can manage to run the code butn just in case I want to ask - is what this function spits out something that can be directly fed into optimize? for instance, do I now replace subdata, submask, and bg in the optimize funcion with trdata, trmask, and trbg?
