"""This package contains weighting functions for computing sound pressure levels 
(SPLs).  The A-weighting function was taken from 
https://gist.github.com/endolith/148112 and used under a BSD license (see the 
preamble to the A-weighting function).  The other weighting functions were 
written in the same form as the A-weighting function.
"""
from __future__ import division
from numpy import pi, polymul
from scipy.signal import bilinear


def A_weighting(fs):
    """Design of an A-weighting filter.
    b, a = A_weighting(fs) designs a digital A-weighting filter for
    sampling frequency `fs`. Usage: y = scipy.signal.lfilter(b, a, x).
    Warning: `fs` should normally be higher than 20 kHz. For example,
    fs = 48000 yields a class 1-compliant filter.
    References:
       [1] IEC/CD 1672: Electroacoustics-Sound Level Meters, Nov. 1996.
    
    Translated from a MATLAB script (which also includes C-weighting,
    octave and one-third-octave digital filters). Author: Christophe
    Couvreur, Faculte Polytechnique de Mons (Belgium)
        couvreur@thor.fpms.ac.be
    Last modification: Aug. 20, 1997, 10:00am. BSD license
    http://www.mathworks.com/matlabcentral/fileexchange/69 Translated
    from adsgn.m to Python 2009-07-14 endolith@gmail.com
    """
    # Definition of analog A-weighting filter according to IEC/CD 1672.
    f1 = 20.598997
    f2 = 107.65265
    f3 = 737.86223
    f4 = 12194.217
    A1000 = 1.9997

    NUMs = [(2*pi * f4)**2 * (10**(A1000/20)), 0, 0, 0, 0]
    DENs = polymul([1, 4*pi * f4, (2*pi * f4)**2],
                   [1, 4*pi * f1, (2*pi * f1)**2])
    DENs = polymul(polymul(DENs, [1, 2*pi * f3]),
                                 [1, 2*pi * f2])

    # Use the bilinear transformation to get the digital filter.
    # (Octave, MATLAB, and PyLab disagree about Fs vs 1/Fs)
    return bilinear(NUMs, DENs, fs)

def M_weighting(fs, flow=40.0, fhigh=300.0):
    """Design of an M-weighting filter for marine bioacoustics
    b, a = M_weighting(fs, flow=40.0, fhigh=300.0) designs a digital M-weighting 
    filter for sampling frequency `fs`. 
    Usage: y = scipy.signal.lfilter(b, a, x).
    Warning: `fs` should normally be higher than 20 kHz.
    Reference:
      Southall, B. L., Bowles, A. E., Ellison, W. T., Finneran, J. J., Gentry, 
      R. L., Greene Jr, C. R., Kastak, D., Ketten, D. R., Miller, J. H., 
      Nachtigall, P. E., Richardson, W. J., Thomas, J. A., and Tyack,
      P. L. (2007), "Marine mammal noise exposure criteria: Initial scientific 
      recommendations," Aquatic Mammals 33(4), 411--521.
    """
    NUMs = [(2*pi * fhigh)**2 + (2*pi * flow)**2, 0, 0]
    DENs = polymul([1, 4*pi * fhigh, (2*pi * fhigh)**2], \
        [1, 4*pi * flow, (2*pi * flow)**2])
    return bilinear(NUMs, DENs, fs)

def weight(data, weighting, fs):
    """Calls the function weighting(fs) to generate a digital filter and then 
    filters data by that filter.  fs is the sampling frequency of the data.
    Usage:
        wdata = weight(data, weighting, fs)
    To provide extra parameters to weighting, use the lambda function. For 
    example:
        wdata = weight(data, lambda fs: M_weighting(fs, flow=150.0, 
                fhigh=1.6e5), fs)
    """
    from scipy.signal import lfilter
    b, a = weighting(fs)
    return(lfilter(b, a, data))
