from __future__ import division
import numpy as np
from scipy.signal import bilinear

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
        
def spl_wav(wavfile, t0=0, t1=-1, weighting='A', tconst=0.125, pref=20.0, \
    cal=1.0, outfile=True, output=False):
    """Return sampled sound pressure levels from data in a WAV file.
    Usage
    
        spl_wav(wavfile, t0=0, t1=-1, weighting='A', tconst=0.125,
            pref=20.0, outfile=True, output=False)
    
    Input Parameters
    
    wavfile: filename of the input WAV file.
    t0: Start time in seconds for data to be analyzed in the WAV file. Defaults 
        to 0.
    t1: End time in seconds for data to be analyzed in the WAV file. A value of 
        -1 (default) represents the end of the file.
    weighting: type of weighting to use. A value of 1 represents no
        weighting. The default value is 'A' for A-weighting. See
        documentation for the function spl().
    tconst: time constant. Defaullts to 0.125 s (fast). This parameter can be 
        the value in seconds or preset values given with the strings 'Fast' 
        (0.125 s), 'Slow' (1.000 s), or "Impulse' (0.035 s).
    pref: reference pressure to normalize output. This should be the numerical 
        value in data is 0 dB. Defaults to 20.0 (micropascal).
    cal: calibration factor of the recording. This is the value that converts 
        data samples to the appropriate pressure units (micro Pa).
    outfile: output file name. A value of True (default) uses the basename of 
        the input file with '_SPL' appended and the extension '.csv'. A string 
        value provides a filename. A value of False results in no output file 
        written.
    output: boolean that determines whether or not the SPL values are returned 
        as output.  A value of True produces output. Defaults to False.
    
    Output
    
    Output is produced only if the parameter output is True.  In this case 
    the output parameters are time, SPL.
        
    time: time values in seconds corresponding SPL samples. The time values are 
        the same as the values in the input WAV file. 
    SPL: a numpy array of sampled sound pressure levels corresponding to the 
        the same sampling times a the elements of data.  Note that the initial 
        elements SPL[i] are based on a truncated history because they only use 
        pressure values from data[i] back to data[0].
    
     The output file is a CSV file with rows of (time, SPL) pairs. 
    """
    from os.path import splitext
    import wavefile
    import csv
    fs, data = wavfile.read(wavfile, t0, t1)
    SPL = spl(data, fs, weighting, tconst, pref, cal)
    time = np.divide(np.arange(SPL.size, dtype=np.float64), np.float64(fs))
    time = np.add(time, np.float64(t0))
    if outfile:
        if outfile == True:
            base, _ = splitext(wavfile)
            outfile = base + '_SPL.csv'
        with open(outfile, 'wb') as csvfile:
            csvwrite = csv.writer(csvfile)
            csvwrite.writerow(['Time (s)', 'SPL (dB)'])
            for n in range(SPL.size):
                csvwrite.writerow([time[n], SPL[n]])
    if output:
        return(time, SPL)

def spl_wav_dir(pattern='*.wav', outname='_SPL.csv', weighting='A.', \
    tconst=0.125, pref=20.0, cal=1.0, icL=False, verbose=False):
    """
    Evaluate all files matching pattern using spl_wav() and write output to 
    files.
    
    Input Parameters
    
    pattern: UNIX-style pathname pattern for input files. See the glob module 
        for details. Defaults to '*.wav', which returns all files in the current 
        directory with extension .wav.
    outname: string appended to the basename of the input file to give the 
        output file name. Defaults to '_SPL.csv'
    weighting: type of weighting to use. A value of 1 represents no
        weighting. The default value is 'A' for A-weighting. See
        documentation for the function spl().
    tconst: time constant. Defaullts to 0.125 s (fast). This parameter can be 
        the value in seconds or preset values given with the strings 'Fast' 
        (0.125 s), 'Slow' (1.000 s), or "Impulse' (0.035 s).
    pref: reference pressure to normalize output. This should be the numerical 
        value in data is 0 dB. Defaults to 20.0 (micropascal).
    cal: calibration factor of the recording. This is the value that converts 
        data samples to the appropriate pressure units (micro Pa).
    icL: boolean that determines if input files were recorded by an icListen 
        hydrophone. If so, the calibration information is extracted from the 
        information embedded in the file. Defaults to False.
    verbose: print progress message. Defaults to False.
    """
    from glob import glob
    from os.path import splitext
    from .recorders.icListen import get_info
    files = glob(pattern)
    fcount = 0
    for file in files:
        if verbose:
            print('Processing file ' + file + ' ...')
        base, _ = splitext(file)
        outfile = base + outname
        if icL:
            cal = get_info(file)['cal']
            if verbose:
                print('Found icListen info. cal: ' + cal)
        spl_wav(file, weighting=weighting, tconst=tconst, pref=pref, \
            cal=cal, outfile=outfile)
        if verbose:
            print('Output written to file ' + outfile)
            fcount += 1
    if verbose:
        print('Complete... Processed ' + str(fcount) + ' files.')

def spl_wav_files(infnames, outname='_SPL.csv', weighting='A', \
    tconst=0.125, pref=20.0, cal=1.0, icL=False, verbose=False):
    """
    Evaluate all files in infnames matching using spl_wav() and write output to 
    files.
    
    Input Parameters
    
    infnames: list of the input files to be processed.
    outname: string appended to the basename of the input file to give the 
        output file name. Defaults to '_SPL.csv'
    weighting: type of weighting to use. A value of 1 represents no
        weighting. The default value is 'A' for A-weighting. See
        documentation for the function spl().
    tconst: time constant. Defaullts to 0.125 s (fast). This parameter can be 
        the value in seconds or preset values given with the strings 'Fast' 
        (0.125 s), 'Slow' (1.000 s), or "Impulse' (0.035 s).
    pref: reference pressure to normalize output. This should be the numerical 
        value in data is 0 dB. Defaults to 20.0 (micropascal).
    cal: calibration factor of the recording. This is the value that converts 
        data samples to the appropriate pressure units (micro Pa).
    icL: boolean that determines if input files were recorded by an icListen 
        hydrophone. If so, the calibration information is extracted from the 
        information embedded in the file. Defaults to False.
    verbose: print progress message. Defaults to False.
    """
    from os.path import splitext
    from .recorders.icListen import get_info
    fcount = 0
    for file in infnames:
        if verbose:
            print('Processing file ' + file + ' ...')
        base, _ = splitext(file)
        outfile = base + outname
        if icL:
            cal = get_info(file)['cal']
            if verbose:
                print('Found icListen info. cal: ' + cal)
        spl_wav(file, weighting=weighting, tconst=tconst, pref=pref, \
            cal=cal, outfile=outfile)
        if verbose:
            print('Output written to file ' + outfile)
            fcount += 1
    if verbose:
        print('Complete... Processed ' + str(fcount) + ' files.')


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

    NUMs = [(2*np.pi * f4)**2 * (10**(A1000/20)), 0, 0, 0, 0]
    DENs = np.polymul([1, 4*np.pi * f4, (2*np.pi * f4)**2],
                   [1, 4*np.pi * f1, (2*np.pi * f1)**2])
    DENs = np.polymul(np.polymul(DENs, [1, 2*np.pi * f3]),
                                 [1, 2*np.pi * f2])

    # Use the bilinear transformation to get the digital filter.
    # (Octave, MATLAB, and PyLab disagree about Fs vs 1/Fs)
    return bilinear(NUMs, DENs, fs)

def M_weighting(fs, flow=40.0, fhigh=300.0):
    """Design of an M-weighting filter for marine bioacoustics 
    b, a = M_weighting(fs, flow=40.0, fhigh=300.0) designs a digital
    M-weighting filter for sampling frequency `fs`. 
    Usage: y = scipy.signal.lfilter(b, a, x).
    Warning: `fs` should normally be higher than 20 kHz.
    Reference:
      Southall, B. L., Bowles, A. E., Ellison, W. T., Finneran, J. J., Gentry, 
      R. L., Greene Jr, C. R., Kastak, D., Ketten, D. R., Miller, J. H., 
      Nachtigall, P. E., Richardson, W. J., Thomas, J. A., and Tyack,
      P. L. (2007), "Marine mammal noise exposure criteria: Initial scientific 
      recommendations," Aquatic Mammals 33(4), 411--521.
    """
    NUMs = [(2*np.pi * fhigh)**2 + (2*np.pi * flow)**2, 0, 0]
    DENs = np.polymul([1, 4*np.pi * fhigh, (2*np.pi * fhigh)**2], \
        [1, 4*np.pi * flow, (2*np.pi * flow)**2])
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
