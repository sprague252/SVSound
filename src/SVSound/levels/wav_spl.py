from __future__ import print_function

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
    from wavesg import wavread
    import numpy as np
    import csv
    fs, data = wavread(wavfile, t0, t1)
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
    from icListen import info
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
    from icListen import info
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
