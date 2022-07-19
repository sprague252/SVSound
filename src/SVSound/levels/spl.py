def spl(data, fs, weighting="A", tconst=0.125, pref=20.0, cal=1.0, pms=0., 
    pms_ret=False):
    """Return an array of sampled sound pressure levels using time constant T.  
    Usage
    
        SPL = spl(data, fs, weighting='A', tconst=0.125, pref=20.0)
    
    Input Parameters
    
    data: an array of sampled sound pressures.
    fs: sampling frequency in hertz.
    weighting: type of weighting to use. This parameter can be a
        string to represent preset values 'A' for A-weighting, 'M'
        for M-weighting (see documentation on the function weight()
        to set frequency parameters). It can also be a function that
        provides digital filter parameters to the weight() function.
        For no weighting, use weighting = 1. The default is 'A'
        weighting.
    tconst: time constant. Defaullts to 0.125 s (fast). This parameter can be 
        the value in seconds or preset values given with the strings 'Fast' 
        (0.125 s), 'Slow' (1.000 s), or "Impulse' (0.035 s).
    pref: reference pressure. Defaults to 20.0 (micropascals, standard for
        atmospheric sounds). Use 1.0 (micropascals) for underwater sounds.
    cal: calibration factor of the recording. This is the value that
        converts data samples to the appropriate pressure units
        (micropascals). The default value is 1 (no calibration
        adjustment).
    pms: an initial value for the mean square pressure 'historical'
        value for time constant. Use this to continue the calculation
        from another recording. Defaults to 0.0.
    pms_return: whether or not to return the mean square pressure value for 
        subsequent calculations. Defaults to False.
    
    Output
    
    SPL: a numpy array of sampled sound pressure levels corresponding to the 
        the same sampling times a the elements of data.  Note that the initial 
        elements SPL[i] are based on a truncated history because they only use 
        pressure values from data[i] back to data[0].
    pms: The mean square sound pressure for use in subsequent
        calculations such as the recording continuing in another
        file. Only returned if the input parameter pms_return is True.
    """
    import numpy as np
    if tconst == 'Fast':
        tconst = 0.125
    elif tconst == 'Slow':
        tconst = 1.0
    elif tconst == 'Impulse':
        tconst = 0.035
    if weighting == 'M':
        weighting = M_weighting
    elif weighting == 'A':
        weighting = A_weighting
    if hasattr(weighting, '__call__'):
        wdata = weight(data, weighting, fs)
    else:
        wdata = data
    caldB = 20. * np.log10(cal)
    spldata = np.zeros(wdata.size)
    for n in range(wdata.size):    
        pms = pms * np.exp(-1 / np.float64(tconst * fs)) + wdata[n]**2 / \
            np.float64(tconst * fs)
        # Calibrate the SPL (better precision with lower numbers).
        spldata[n] = (10.0 * np.log10(pms/pref**2)  + caldB)
    if pms_ret:
        return(spldata, pms)
    else:
        return(spldata)
