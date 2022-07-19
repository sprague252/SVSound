from __future__ import division
    
def sel(data, fs, weighting="A", pref=20.0, cal=1.0, amb=-1.0e3, pms=0., \
    pms_ret=False):
    """Return the sound exposure level of a sampled sound.  
    Usage
    
        SEL = sel(data, fs, weighting='A', pref=20.0, cal=1.0,
                  amb=-1.0e3, pms=0., pms_ret=False))
    
    Input Parameters
    
    data: an array of sampled sound pressures.
    fs: sampling frequency in hertz.
    weighting: type of weighting to use. This parameter can be a
        string to represent preset values 'A' for A-weighting, 'M'
        for M-weighting (see documentation on the function
        SPLWeighting.weight() to set frequency parameters). The
        parameter weighting can also be a function that provides
        digital filter parameters to the SPLWeighting.weight()
        function. For no weighting, use weighting=1. The default is
        'A' weighting.
    pref: reference pressure. Defaults to 20.0 (micropascals, the
        default value for atmospheric sounds). Use 1.0 (micropascals)
        for underwater sounds.
    cal: calibration factor of the recording. This is the value that converts 
        data samples to the appropriate pressure units (micro Pa).
    pms: an initial value for the mean square pressure 'historical' value for 
        time constant. Defaults to 0.0.
    pms_return: whether or not to return the mean square pressure
        value (in linear units) for subsequent calculations.  This
        allows the SEL to be calculated from several sound
        recordings. Defaults to False.
    
    Output
    
    SEL: the sound exposure level of the recording in dB re pref.
        Note that if pms_return is True, the intermediate value of
        mean square pressure is returned in linear units (see
        pms_return above). 
    """
    import numpy as np
    from . import SPLWeighting
    p2amb = pref**2 * 10.0**(amb/10.0)
    if weighting == 'M':
        weighting = SPLWeighting.M_weighting
    elif weighting == 'A':
        weighting = SPLWeighting.A_weighting
    if hasattr(weighting, '__call__'):
        wdata = SPLWeighting.weight(data, weighting, fs)
    else:
        wdata = data
    for p in wdata:
        pms = pms + ((p * cal)**2 - p2amb) / fs
    if pms_ret:
        return(pms)
    else:
        selcal = 10.0 * np.log10(pms / pref**2)
        return(selcal)
