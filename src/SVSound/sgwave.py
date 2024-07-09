"""
Sgwave is a package that processes and stores spectrogram data in HDF5
format data files. There are functions to process WAV files and store
the spectrograms in the data files and functions to extract spectrograms
from the data files to perform various analyses.

2017-09-26 23:05
"""
from __future__ import division

def sgwavfile(wavfile, t0=0, t1=-1, nfft=2048, noverlap=-1, window='hann', \
    flength=60., ftime=True, datestring=None, ftype=None, verbose=False):
    """Compute spectrograms of a .WAV file and save it to HDF5 files
    (with bzip2 compression implemented). This function uses
    scipy.signal.spectrogram to compute the spectrograms of the data.
    
    sgwavfile(wavfile, t0, t1, nfft, noverlap, window, flength, ftime, ftype, verbose)
    
    Input Parameters
    wavfile: filename of the WAV file to process
    t0: start time (in seconds) for processing. The default value is 0.
    t1: end time (in seconds) for processing. Use -1 for the end of
        the file. The default value is -1.
    nfft: number of points in each FFT. The default value is 2048.
    noverlap: the number of waveform points to overlap in each FFT.
        Use -1 for nfft/2. The default value is -1.
    window: FFT window function or function name as required by
        scipy.signal.spectrogram. The default value is 'hann' (for a
        Hanning window).
    flength: length of each spectrogram file. May either be in time
        (seconds) or samples. The default value is 60 seconds.
    ftime: boolean that specifies whether flength is a time (True) or
        a number of samples (False). The default value is True.
    datestring: set the datestring for the sprctrogram file. If
        datestring is None the datestring is extracted from the WAV
        file. The default value is None.
    ftype: The type of WAV file. This parameter allows the extraction
        of proprietary metadata recorded in some WAV files by
        different recorder. Use None for a standard WAV file,
        'icListen' for an icListen file, and 'zoom' for a Zoom
        recorder WAV file. The default value is False.
    verbose: enable verbose output. The default value is None.
    
    Output
    The output files use the same name as the input files with a
    "_n_sg.h5." with (n=file number) appended before the file extension.
    If the input file is 'fname.wav', and there are 10 spectrogram files
    generated from it, the output files will be named 'fname_0_sg.h5'
    through 'fname_9_sg.h5'. The sgwavfile function returns a value of 0
    upon successful completion.
    """
    import numpy as np
    import scipy.signal as sig
    import tables
    from os.path import exists, basename, dirname
    from os import stat
    from SVSound import wavefile
    import logging
    ch = logging.StreamHandler()
    formatter = logging.Formatter('%(asctime)s %(levelname)8s %(name)s | %(message)s')
    ch.setFormatter(formatter)
    logger = logging.getLogger('sgwavfile')
    logger.addHandler(ch)
    if verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)
    if (noverlap == -1):
        noverlap = nfft//2
    # Store each sample at a float64.
    atom = tables.Float64Atom()
    # This is the compression algorithm for the output file.
    filters = tables.Filters(complevel=9, complib='bzip2', shuffle=True, \
                fletcher32=True)
    with open(wavfile, 'rb') as infile:
        info = wavefile.getInfo(infile, ftype)
        if datestring == None:
            # Get the datestring from the WAV file.
            if (info['wavetype'] == 'zoom'):
                # Time information is stored in the WAV file in the time
                # zone for which the recorder was set.
                datestring = (info['OriginationDate'] + 'T' +
                    info['OriginationTime'])
            elif (info['wavetype'] == 'icListen'):
                # Extract the time and date information from the
                # filename.
                tdinfo = wavfile.replace('.', '_').split('_')[1:3]
                datestring = (tdinfo[0][:4]+ '-' + tdinfo[0][4:6] +
                    '-' + tdinfo[0][6:8] + 'T' + tdinfo[1][:2] + ':'
                    + tdinfo[1][2:4] + ':' + tdinfo[1][4:6] +
                    'Z')
            elif (info['wavetype'] == 'decimus'):
                # Extract the time and date information from the
                # filename.
                tdinfo = wavfile.replace('.', '_').split('_')[1:4]
                datestring = (tdinfo[0][:4]+ '-' + tdinfo[0][4:6] +
                    '-' + tdinfo[0][6:8] + 'T' + tdinfo[1][:2] + ':'
                    + tdinfo[1][2:4] + ':' + tdinfo[1][4:6] + '.' +
                    tdinfo[2] + 'Z')
            else:
                # Try to get datestring from WAV file modification time.
                fm = stat(wavfile).st_mtime
                datestring = (str( 
                    np.datetime64(np.int64(fm), 's') - 
                        np.timedelta64( 
                            np.int64( 
                                np.round( 
                                    info['Nsamples'] / info['fs'] 
                                ) 
                            ), 
                            's' 
                        ) 
                    ) + 
                    'Z')
        # flength determines the size of the input file slice in seconds
        # (ftime==True) or samples (ftime ==False) for each output file.
        # Note that there will also be some overlap samples at the end
        # of intermediate input files.
        if ftime==False:
            # flength is in samples. Convert to seconds for
            # compatibility with other code.
            flength = flength / info['fs']
        # n0 is the initial sample to process.
        n0 = int(t0 * info['fs'])
        # n1 is the last sample to process.
        if (t1 == -1):
            n1 = info['Nsamples'] - 1
        else:
            n1 = int(t1 * info['fs']) - 1
        # Determine the number of samples to process.
        nn = n1 - n0 + 1
        # Determine the number of output files. There may be extra
        # samples at the end.
        nfiles = np.int64(np.floor(nn / info['fs'] / flength))
        # nnfile is the number of wave samples in each file, not
        # including overlap.
        nnfile = np.int64(flength * info['fs'])
        # nskip is the number of samples to skip for each FFT.
        nskip = nfft - noverlap
        # nnfft is the number of FFTs to process in each complete input
        # file.
        nnfft = np.int64(np.ceil((flength * info['fs'] - nfft) /
            nskip + 1))
        # nend is the number of samples from the first input file sample
        # to the last input file sample to process in each output file.
        # There could be some extra samples past flength to complete the
        # last FFT.
        nend = (nnfft - 1) * nskip + nfft - 1
        wavdir = dirname(wavfile)
        wavbase = basename(wavfile)
        if (wavdir == ''):
            outdir = ''
        else:
            outdir = wavdir + '/'
        for nfile in range(nfiles - 1):
            # Name of output file.
            outfile = outdir + wavbase.rsplit('.', 1)[0] + \
                '_{0:=-02d}_sg.h5'.format(np.int64(t0 / flength) + \
                nfile)
            logger.info('Infile: {}, Outfile: {}'.format(wavfile,
                outfile))
            with tables.open_file(outfile, mode='w', filters=filters) \
                as fout:
                out_t = fout.create_earray(fout.root, 't', atom,
                    (0,), title='Time (s)', expectedrows=nnfft)
                out_f = fout.create_earray(fout.root, 'f', atom,
                    (0,), title='Frequency (Hz)', expectedrows=nnfft)
                out_sg = fout.create_earray(fout.root, 'sg', atom,
                    (nfft//2+1, 0), title='Spectrogram',
                    expectedrows=nnfft)
                out_sg.attrs.inputfile = wavfile
                out_t.attrs.datestring = datestring
                out_sg.attrs.fs = info['fs']
                out_sg.attrs.nfft = nfft
                out_sg.attrs.noverlap = noverlap
                out_sg.attrs.window = window
                # The following is right because nend is the number of
                # input file samples used, not an absolute sample
                # number.
                out_sg.attrs.Nsamples = nend + 1
                t0f = t0 + (nfile * nnfile) / info['fs']
                logger.info('Infile: {}, t0f: {}'.format(wavfile,
                    t0f))
                t1f = t0f + nend / info['fs']
                logger.info('Infile: {}, t1f: '.format(wavfile, t1f))
                wave = wavefile.wave_chunk(infile, info, t0 = t0f, 
                    t1 = t1f, verbose=verbose)
                logger.info(('Infile: {}, Returned from '
                    + 'wave_chunk').format(wavfile)) 
                # Apply the calibration value if it exists.
                if ('cal' in info):              
                    f, t, sg = sig.spectrogram(np.multiply(wave,
                        info['cal']), info['fs'], window='hann',
                        nperseg=nfft, noverlap=noverlap)
                # Else, use the waveform as is in file. 
                else:
                    f, t, sg = sig.spectrogram(wave, info['fs'],
                        window='hann', nperseg=nfft,
                        noverlap=noverlap)                
                out_t.append(t + t0f)
                out_t.flush()
                out_f.append(f)
                out_f.flush()
                out_sg.append(sg)
                out_sg.flush()
        # Take care of last file so it does not pass the end of the
        # input file.
        nfile = nfiles - 1
        outfile = outdir + wavbase.rsplit('.', 1)[0] + \
            '_{0:=-02d}_sg.h5'.format(np.int64(t0 / flength) + \
            nfile)
        logger.info('Infile: {}, Outfile: {}'.format(wavfile,
            outfile))
        with tables.open_file(outfile, mode='w', filters=filters) as \
            fout:
            out_t = fout.create_earray(fout.root, 't', atom, (0,),
                title='Time (s)', expectedrows=nnfft)
            out_f = fout.create_earray(fout.root, 'f', atom, (0,), 
                title='Frequency (Hz)', expectedrows=nnfft)
            out_sg = fout.create_earray(fout.root, 'sg', atom,
                (nfft//2+1, 0), title='Spectrogram',
                expectedrows=nnfft)
            out_sg.attrs.inputfile = wavfile
            out_t.attrs.datestring = datestring
            out_sg.attrs.fs = info['fs']
            out_sg.attrs.nfft = nfft
            out_sg.attrs.noverlap = noverlap
            out_sg.attrs.window = window
            out_sg.attrs.Nsamples = nend + 1
            t0f = t0 + (nfile * nnfile) / info['fs']
            logger.info('Infile: {}, t0f: {}'.format(wavfile, t0f))
            t1f = t0f + nend / info['fs']
            # If the extra samples to complete the last FFT exceed the
            # end of the file, reduce neng by nskip and adjust.
            if (t1f * info['fs'] > info['Nsamples']):
                nend = nend - nskip
                t1f = t0f + nend / info['fs']
            logger.info('Infile: {}, t1f: '.format(wavfile, t1f))
            wave = wavefile.wave_chunk(infile, info, t0=t0f, t1=t1f,
                verbose=verbose)
            logger.info(('Infile: {}, Returned from '
                + 'wave_chunk').format(wavfile)) 
            # Apply the calibration value if it exists.
            if ('cal' in info):              
                f, t, sg = sig.spectrogram(np.multiply(wave,
                    info['cal']), info['fs'], window='hann',
                    nperseg=nfft, noverlap=noverlap)
            # Else, use the waveform as is in file. 
            else:
                f, t, sg = sig.spectrogram(wave, info['fs'],
                    window='hann', nperseg=nfft,
                    noverlap=noverlap)                
            out_t.append(t + t0f)
            out_t.flush()
            out_f.append(f)
            out_f.flush()
            out_sg.append(sg)
            out_sg.flush()
        # Take care of samples at the end.
        nnlast = nn - nfiles * nnfile
        if (nnlast >= nfft):
            # nnfft is the number of FFTs to process in the input file.  
            # Using integer math to truncate.
            nnfftlast = np.int64(np.floor((nnlast - nfft) / nskip + 1))
            # nend is the last input file sample to process.
            nend = (nnfftlast - 1) * nskip + nfft - 1
            outfile = wavfile.rsplit('.', 1)[0] + \
                '_{0:=-02d}_sg.h5'.format(np.int64(t0 / flength) + \
                nfiles)
            logger.info('Infile: {}, Outfile: {}'.format(wavfile,
                outfile))
            with tables.open_file(outfile, mode='w', filters=filters) \
                as fout:
                out_t = fout.create_earray(fout.root, 't', atom,
                    (0,), title='Time (s)', expectedrows=nnfftlast)
                out_f = fout.create_earray(fout.root, 'f', atom,
                    (0,), title='Frequency (Hz)',
                    expectedrows=nnfftlast)
                out_sg = fout.create_earray(fout.root, 'sg', atom,
                    (nfft//2+1, 0), title='Spectrogram',
                    expectedrows=nnfftlast)
                out_t.attrs.datestring = datestring
                out_t.attrs.unit = 's'
                out_f.attrs.unit = 'Hz'
                out_sg.attrs.inputfile = wavfile
                out_sg.attrs.fs = info['fs']
                out_sg.attrs.nfft = nfft
                out_sg.attrs.noverlap = noverlap
                out_sg.attrs.window = window
                out_sg.attrs.Nsamples = nend + 1
                out_sg.attrs.unit = chr(956) + 'Pa^2/Hz'
                out_sg.attrs.composite = False
                t0f = t0 + ((nfiles - 1) * nnfile) / info['fs']
                logger.info('Infile: {}, t0f: {}'.format(wavfile,
                    t0f))
                t1f = t0f + nend / info['fs']
                logger.info('Infile: {}, t1f: {}'.format(wavfile,
                    t1f))
                wave = wavefile.wave_chunk(infile, info, t0=t0f,
                    t1=t1f, verbose=verbose)
#               logger.info('wave read. wave.size: ', wave.size)
                # Apply the calibration value if it exists.
                if ('cal' in info):              
                    f, t, sg = sig.spectrogram(np.multiply(wave,
                        info['cal']), info['fs'], window='hann',
                        nperseg=nfft, noverlap=noverlap)
                # Else, use the waveform as is in file. 
                else:
                    f, t, sg = sig.spectrogram(wave, info['fs'],
                        window='hann', nperseg=nfft,
                        noverlap=noverlap)                
#               logger.info('sg complete.')
#               logger.info('  f.size: ', f.size)
#               logger.info('  t.size: ', t.size)
#               logger.info('  sg.size: ', sg.size)
#               logger.info('')
                out_t.append(t + t0f)
                out_t.flush()
                out_f.append(f)
                out_f.flush()
                out_sg.append(sg)
                out_sg.flush()
            logger.info('Infile: {}, Processing ' + \
                'complete.'.format(wavfile))
            logger.info(('Infile: {}, {:d} input file samples ' + \
                'processed.').format(wavfile, nn))
            logger.info(('Infile: {}, {:d} files ' 
                + 'written').format(wavfile, nfiles + 1))
            logger.info(('Infile: {}, {:d} FFTs in files 0 -'
                + ' {:d}').format(wavfile, nnfft, nfiles - 1))
            logger.info(('Infile: {}, {:d} FFTs in file '
                + '{:d}').format(wavfile, nnfftlast,                    
                nfiles))
        else:
            logger.info(('Infile: {}, Processing '
                + 'complete.').format(wavfile))
            logger.info(('Infile: {}, {:d} input file samples '
                + 'processed.').format(wavfile, nn))
            logger.info(('Infile: {}, {:d} files '
                + 'written').format(wavfile, nfiles))
            logger.info(('Infile: {}, {:d} FFTs in each '
                + 'file').format(wavfile, nnfft))
    # Return a value of 0 for successful completion.
    logger.removeHandler(ch)
    return(0)

def get_sginfo(sgfile):
    import numpy as np
    import tables
    with tables.open_file(sgfile, mode='r') as infile:
        sginfo = {}
        if ('fs' in infile.root.sg.attrs._f_list()):
            sginfo['fs'] = infile.root.sg.attrs.fs
        if ('nfft' in infile.root.sg.attrs._f_list()):
            sginfo['nfft'] = infile.root.sg.attrs.nfft
        if ('fs' in infile.root.sg.attrs._f_list()):
            sginfo['noverlap'] = infile.root.sg.attrs.noverlap
        if ('window' in infile.root.sg.attrs._f_list()):
            sginfo['window'] = infile.root.sg.attrs.window
        if ('datestring' in infile.root.t.attrs._f_list()):
            sginfo['datestring'] = infile.root.t.attrs.datestring
        if ('unit' in infile.root.t.attrs._f_list()):
            sginfo['t unit'] = infile.root.t.attrs.unit
        if ('unit' in infile.root.f.attrs._f_list()):
            sginfo['f unit'] = infile.root.f.attrs.unit
        if ('unit' in infile.root.sg.attrs._f_list()):
            sginfo['sg unit'] = infile.root.sg.attrs.unit
        if ('WAV file' in infile.root.sg.attrs._f_list()):
            sginfo['WAV file'] = infile.root.sg.attrs.inputfile
        if ('WAV nsamp' in infile.root.sg.attrs._f_list()):
            sginfo['WAV nsamp'] = infile.root.sg.attrs.Nsamples
        if ('composite' in infile.root.sg.attrs._f_list()):
            sginfo['composite'] = infile.root.sg.attrs.composite
        else:
            sginfo['composite'] = False
    return(sginfo)

def sgread(sgfile, t0=0, t1=-1, f0=0, f1=-1):
    """Read an HDF5 file written by sgwavefile and extract the spectrogram information \
    from it.
    """
    import numpy as np
    import tables
    with tables.open_file(sgfile, mode='r') as infile:
        t = infile.root.t.read()
        f = infile.root.f.read()
        if (t.size == 0):
            return([], [], [])
        if (t1 == -1):
            t1 = t[-1]
        if (f1 == -1):
            f1 = f[-1]
        tb = np.logical_and(t >= t0, t <= t1)
        fb = np.logical_and(f >= f0, f <= f1)
        tout = t[tb]
        t_ind = np.array(np.where(tb))
        t_ind0 = t_ind.min()
        t_ind1 = t_ind.max()
        fout = f[fb]
        f_ind = np.array(np.where(fb))
        f_ind0 = f_ind.min()
        f_ind1 = f_ind.max()
        sgout = infile.root.sg[t_ind0 : t_ind1 + 1, f_ind0 : f_ind1+1]
    return(fout, tout, sgout)

def sgreadsample(sgfile, nt0=0, nt1=-1, nf0=0, nf1=-1):
    """Read an HDF5 file written by sgwavefile and extract the spectrogram information \
    from it on a given time and frequency sample range.
    """
    import numpy as np
    import tables
    with tables.open_file(sgfile, mode='r') as infile:
        t = infile.root.t.read()
        f = infile.root.f.read()
        if (t.size == 0):
            return(np.array([]), np.array([]), np.array([]))
        if (nt1 == -1):
            nt1 = t.size
        if (nf1 == -1):
            nf1 = f.size
        if nt0 > t.size:
            try:
                raise ValueError('Index nt0 > maximum time sample')
            except Exception as e:
                e.add_note('nt0: {:}'.format(nt0))
                e.add_note('t.size: {:}'.format(t.size))
                raise
        if nt1 > t.size:
            try:
                raise ValueError('Index nt1 > maximum time sample')
            except Exception as e:
                e.add_note('nt1: {:}'.format(nt1))
                e.add_note('t.size: {:}'.format(t.size))
                raise
        if nf0 > f.size:
            try:
                raise ValueError('Index nf0 > maximum frequency sample')
            except Exception as e:
                e.add_note('nf0: {:}'.format(nf0))
                e.add_note('f.size: {:}'.format(f.size))
                raise
        if nf1 > f.size:
            try:
                raise ValueError('Index nf1 > maximum frequency sample')
            except Exception as e:
                e.add_note('nf1: {:}'.format(nf1))
                e.add_note('f.size: {:}'.format(f.size))
                raise
        fout = f[nf0:nf1]
        tout = t[nt0:nt1]      
        sgout = infile.root.sg[nt0:nt1, nf0:nf1]
    return(fout, tout, sgout)

    
def compsg(sgfiles, t0=0, t1=-1, avtime=-1, savefile=False, foutname=None):
    """Produce a composite spectrogram from data files written by sgwavfile.
    """
    import numpy as np
    from os.path import exists, basename, dirname
    if isinstance(sgfiles, str):
        # There is only one file not in a list. Put it in a list.
        sgfiles = [sgfiles]
    elif not (isinstance(sgfiles, list) or isinstance(sgfiles, tuple)):
        raise ValueError('Parameter not list, tuple, or string', sgfiles)
    taps = np.array([], dtype=np.float64)
    datet = np.array([], dtype='M8[ns]')
    # Process first file.
    info = get_sginfo(sgfiles[0])
    if (type(info['datestring']) == np.bytes_) or \
        (type(info['datestring']) == bytes):
        info['datestring'] = info['datestring'].decode('utf-8')
    # Note: We are assuming that the utc value for all files is the same.
    if (info['datestring'][-1] == 'Z'):
        # datestring is in UTC. Avoid depreciation warning in np.datetime64.
        utc = True
        datet0 = np.datetime64(info['datestring'][:-1])
    else:
        # datestring is in local time.
        utc = False
        datet0 = np.datetime64(info['datestring'])
    nskip = info['nfft'] - info['noverlap']
    dt = nskip / info['fs']
    if (len(sgfiles) == 1):
        f, t, sg = sgread(sgfiles[0], t0=t0, t1=t1)
    else:
        f, t, sg = sgread(sgfiles[0], t0=t0, t1=-1)
    t0f = t[0] - dt
    # Determine the number of average power spectra in for the file.
    if (avtime == -1):
        taps = np.append(taps, (t[0] + t[-1]) / 2)
        datet = np.append(datet, datet0 + \
            np.timedelta64(np.int64(taps[-1] * 1.e6), 'us'))
        psav = np.sum(sg, axis=1).reshape(f.size, 1) / t.size
        csg = psav
    else:
        # Find the number of spectrograms for which the last spectrogram is based on
        # samples entirely within the avtime window.
        avnn = np.int64(np.floor((avtime * info['fs'] - info['nfft']) / 
            (info['nfft'] - info['noverlap']) + 1))
        # nnpsav is the number of average power spectra in the file.
        nnpsav = np.int64(np.floor(t.size / avnn))
        for n in range(nnpsav):
            nst0 = n * avnn
            nst1 = (n + 1) * avnn
            #if (nst1 >= t.size):
            #    nst1 = t.size
            taps = np.append(taps, (t[nst0] + t[nst1]) / 2.0)
            datet = np.append(datet, datet0 + \
                np.timedelta64(np.int64(taps[-1] * 1.e6), 'us'))
            psav = np.sum(sg[:, nst0:nst1], axis=1).reshape(f.size, 1) / (nst1 - nst0)
            if (n == 0):
                csg = psav
            else:
                csg = np.append(csg, psav, axis=1)
    # Process the intermediate files.
    for sgfile in sgfiles[1:-1]:
        info = get_sginfo(sgfile)
        if (type(info['datestring']) == np.bytes_) or \
            (type(info['datestring']) == bytes):
            info['datestring'] = info['datestring'].decode('utf-8')

        if (info['datestring'][-1] == 'Z'):
            # datestring is in UTC. Avoid depreciation warning in np.datetime64.
            datet0 = np.datetime64(info['datestring'][:-1])
        else:
            # datestring is in local time.
            datet0 = np.datetime64(info['datestring'])
        nskip = info['nfft'] - info['noverlap']
        dt = nskip / info['fs']
        f, t, sg = sgread(sgfile, t0=0, t1=-1)
        t0f = t[0] - dt
        # Determine the number of average power spectra in for the file.
        if (avtime == -1):
            taps = np.append(taps, (t[0] + t[-1]) / 2.0)
            datet = np.append(datet, datet0 + \
                np.timedelta64(np.int64(taps[-1] * 1.e6), 'us'))
            psav = np.sum(sg, axis=1).reshape(f.size, 1) / t.size
            csg = np.append(csg, psav, axis=1)
        else:
            # Find the number of spectrograms for which the last spectrogram is based on
            # samples entirely within the avtime window.
            avnn = np.int64(np.floor((avtime * info['fs'] - info['nfft']) / 
                (info['nfft'] - info['noverlap']) + 1))
            nnpsav = np.int64(np.floor(t.size / avnn))
            for n in range(nnpsav):
                nst0 = n * avnn
                nst1 = (n + 1) * avnn
                #if (nst1 >= t.size):
                #    nst1 = t.size
                taps = np.append(taps, (t[nst0] + t[nst1]) / 2.0)
                datet = np.append(datet, datet0 + \
                    np.timedelta64(np.int64(taps[-1] * 1.e6), 'us'))
                psav = np.sum(sg[:, nst0:nst1], axis=1).reshape(f.size, 1) / \
                    (nst1 - nst0)
                csg = np.append(csg, psav, axis=1)
    if (len(sgfiles) > 1):
        # Process the last file.
        info = get_sginfo(sgfiles[-1])
        if (type(info['datestring']) == np.bytes_) or \
            (type(info['datestring']) == bytes):
            info['datestring'] = info['datestring'].decode('utf-8')
        if (info['datestring'][-1] == 'Z'):
            # datestring is in UTC. Avoid depreciation warning in np.datetime64.
            datet0 = np.datetime64(info['datestring'][:-1])
        else:
            # datestring is in local time.
            datet0 = np.datetime64(info['datestring'])
        nskip = info['nfft'] - info['noverlap']
        dt = nskip / info['fs']
        f, t, sg = sgread(sgfiles[-1], t0=0, t1=t1)
        t0f = t[0] - dt
        # Determine the number of average power spectra in for the file.
        if (avtime == -1):
            taps = np.append(taps, (t[0] + t[-1]) / 2.0)
            datet = np.append(datet, datet0 + \
                np.timedelta64(np.int64(taps[-1] * 1.e6), 'us'))
            psav = np.sum(sg, axis=1).reshape(f.size, 1) / t.size
            csg = np.append(csg, psav, axis=1)
        else:
            # Find the number of spectrograms for which the last spectrogram is based on
            # samples entirely within the avtime window.
            avnn = np.int64(np.floor((avtime * info['fs'] - info['nfft']) / 
                (info['nfft'] - info['noverlap']) + 1))
            # nnpsav is the number of average power spectra in the file.
            nnpsav = np.int64(np.floor(t.size / avnn))
            for n in range(nnpsav):
                nst0 = n * avnn
                nst1 = (n + 1) * avnn
                #if (nst1 >= t.size):
                #    nst1 = t.size
                datet = np.append(datet, datet0 + \
                    np.timedelta64(np.int64(taps[-1] * 1.e6), 'us'))
                taps = np.append(taps, (t[nst0] + t[nst1]) / 2.0)
                psav = np.sum(sg[:, nst0:nst1], axis=1).reshape(f.size, 1) / \
                    (nst1 - nst0)
                csg = np.append(csg, psav, axis=1)
    if savefile:
        if (not foutname):
            foutname = basename(sgfiles[0]).rsplit('_sg', 1)[0] + \
                '_csg.h5'
        compsg_save(foutname, csg, f, taps, datet, utc, info)
    return(f, taps, datet, utc, csg)

def compsg_save(fname, csg, f, t, datet, utc=True, sginfo=False):
    import numpy as np
    import tables
    filters = tables.Filters(complevel=9, complib='bzip2', shuffle=True, \
                fletcher32=True)
    atom = tables.Float64Atom()
    with tables.open_file(fname, mode='w', filters=filters) as fout:
        out_t = fout.create_earray(fout.root, 't', atom, (0,), \
            title='Time (s)', expectedrows=t.size)
        out_datet = fout.create_earray(fout.root, 'datet', atom, (0,), \
            title='Time (s)', expectedrows=datet.size)
        out_f = fout.create_earray(fout.root, 'f', atom, (0,), \
            title='Frequency (Hz)', expectedrows=f.size)        
        out_sg = fout.create_earray(fout.root, 'sg', atom, (csg.shape[0], 0), \
            title='Spectrogram', expectedrows=csg.shape[1])
        if (sginfo):
            sgkeys = sginfo.keys()            
            if 'fs' in sgkeys:
                out_sg.attrs.fs = sginfo['fs']
            if 'nfft' in sgkeys:
                 out_sg.attrs.nfft = sginfo['nfft']
            if 'noverlap' in sgkeys:
                out_sg.attrs.noverlap = sginfo['noverlap']
            if 'window' in sgkeys:
                out_sg.attrs.window = sginfo['window']
            if 't unit' in sgkeys:
                out_t.attrs.unit = sginfo['t unit']
            if 'f unit' in sgkeys:
                out_f.attrs.unit = sginfo['f unit']
            if 'sg unit' in sgkeys:
                out_sg.attrs.unit = sginfo['sg unit']
        out_t.attrs.datestring = str(datet[0])
        out_sg.attrs.composite = True
        out_t.append(t)
        out_t.flush()
        out_datet.append(datet)
        out_datet.flush()
        out_f.append(f)
        out_f.flush()
        out_sg.append(csg)
        out_sg.flush()
        
def compsg_plot_files(sgfiles, t0=0, t1=-1, avtime=-1, flim=False, dtformat=None,
    savefile=False, foutname=None):
    """Use compsg to produce a composite spectrogram and plot it.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from datetime import datetime
    from os.path import basename, dirname
    import warnings
    f, _, datet, utc, csg = compsg(sgfiles, t0, t1, avtime)
    # Need to convert from np.datetime64 to datetime.datetime to plot
    # with matplotlib.
    dt0 = _dt_convert(datet[0])
    dt1 = _dt_convert(datet[-1])
    x_lims = mdates.date2num([dt0, dt1])
    fig, ax = plt.subplots()
    sg = plt.imshow(10*np.log10(csg), origin='lower', aspect='auto',
        extent=[x_lims[0], x_lims[-1], f[0], f[-1]])
    ax.xaxis_date()
    if (dtformat):
        date_format = mdates.DateFormatter(dtformat)
    else:
        date_format = mdates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(date_format)
    if (utc):
        ax.set_xlabel('UTC Time (HH:MM)')
    else:
        ax.set_xlabel('Time (HH:MM)')
    cb = plt.colorbar(shrink=0.5)
    cb.set_label('Power Spectral Density (dB)')
    if (flim):
      ax.set_ylim(flim)
    fig.tight_layout()
    if (savefile):
        if (not foutname):
            foutname = basename(sgfiles[0]).rsplit('_sg', 1)[0] + \
                '_csg.pdf'
        fig.savefig(foutname)
        plt.close(fig)
    else:
        fig.show()
        return(fig, ax, sg, cb)
    

def psbsum(sgfiles, f0=0, f1=-1, t0=0, t1=-1, avtime=-1,
    savefile=False, foutname=None):
    """Produce a power spectral band sum from data files written by
    sgwavfile.
    """
    import numpy as np
    if isinstance(sgfiles, str):
        # There is only one file not in a list. Put it in a list.
        sgfiles = [sgfiles]
    elif not (isinstance(sgfiles, list) or isinstance(sgfiles, tuple)):
        raise ValueError('Parameter not list, tuple, or string', sgfiles)
    taps = np.array([], dtype=np.float64)
    datet = np.array([], dtype='M8[ns]')
    psbsum = np.array([], dtype=np.float64)
    # Process first file.
    info = get_sginfo(sgfiles[0])
    if (type(info['datestring']) == np.bytes_) or \
        (type(info['datestring']) == bytes):
        info['datestring'] = info['datestring'].decode('utf-8')
    # Note: We are assuming that the utc value for all files is the same.
    if (info['datestring'][-1] == 'Z'):
        # datestring is in UTC. Avoid depreciation warning in np.datetime64.
        utc = True
        datet0 = np.datetime64(info['datestring'][:-1])
    else:
        # datestring is in local time.
        utc = False
        datet0 = np.datetime64(info['datestring'])
    nskip = info['nfft'] - info['noverlap']
    dt = nskip / info['fs']
    if (len(sgfiles) == 1):
        f, t, sg = sgread(sgfiles[0], t0=t0, t1=t1)
    else:
        f, t, sg = sgread(sgfiles[0], t0=t0, t1=-1)
    t0f = t[0] - dt
    # nf0 is the index of the lowest frequency band to be included in the sum.
    nf0 = np.int64(np.floor(f0 * info['nfft'] / info['fs']))
    # nf1 is the index of the highest frequency band to be included in the sum.
    if (f1==-1):
        f1 = f[-1]
    nf1 = np.int64(np.floor(f1 * info['nfft'] / info['fs']))
    # Determine the number of average power spectra in for the file.
    if (avtime == -1):
        taps = np.append(taps, (t[0] + t[-1]) / 2.0)
        datet = np.append(datet, datet0 + \
            np.timedelta64(np.int64(taps[-1] * 1.e6), 'us'))
        psbsum = np.append(psbsum, np.sum(sg[nf0:nf1+1].flatten()) / t.size)
    else:
        # Find the number of spectrograms for which the last spectrogram is based on
        # samples entirely within the avtime window.
        avnn = np.int64(np.floor((avtime * info['fs'] - info['nfft']) / 
            nskip + 1))
        # nnpsav is the number of average power spectra in the file.
        nnpsav = t.size // avnn
        for n in range(nnpsav):
            nst0 = n * avnn
            nst1 = (n + 1) * avnn
            if (nst1 >= t.size):
                nst1 = t.size
            taps = np.append(taps, (t[nst0] + t[nst1-1]) / 2.0)
            datet = np.append(datet, datet0 + \
                np.timedelta64(np.int64(taps[-1] * 1.e6), 'us'))
            psbsum = np.append(psbsum, np.sum(sg[nf0:nf1+1, \
                nst0:nst1].flatten()) / (nst1 - nst0))
    # Process the intermediate files.
    for sgfile in sgfiles[1:-1]:
        info = get_sginfo(sgfile)
        if (type(info['datestring']) == np.bytes_) or \
            (type(info['datestring']) == bytes):
            info['datestring'] = info['datestring'].decode('utf-8')
        if (info['datestring'][-1] == 'Z'):
            if (not utc):
                warnings.warn('Inconsitent UTC time encoding') 
            # datestring is in UTC. Avoid depreciation warning in np.datetime64.
            datet0 = np.datetime64(info['datestring'][:-1])
        else:
            if (utc):
                warnings.warn('Inconsitent UTC time encoding') 
            # datestring is in local time.
            datet0 = np.datetime64(info['datestring'])
        nskip = info['nfft'] - info['noverlap']
        dt = nskip / info['fs']
        f, t, sg = sgread(sgfile, t0=0, t1=-1)
        t0f = t[0] - dt
        # nf0 is the index of the lowest frequency band to be included in the sum.
        nf0 = np.int64(np.floor(f0 * info['nfft'] / info['fs']))
        # nf1 is the index of the highest frequency band to be included in the sum.
        nf1 = np.int64(np.floor(f1 * info['nfft'] / info['fs']))
        # Determine the number of average power spectra in for the file.
        if (avtime == -1):
            taps = np.append(taps, (t[0] + t[-1]) / 2.0)
            datet = np.append(datet, datet0 + \
                np.timedelta64(np.int64(taps[-1] * 1.e6), 'us'))
            psbsum = np.append(psbsum, \
                np.sum(sg[nf0:nf1+1].flatten()) / t.size)
        else:
            # Find the number of spectrograms for which the last spectrogram is based on
            # samples entirely within the avtime window.
            avnn = np.int64(np.floor((avtime * info['fs'] - info['nfft']) / 
                nskip + 1))
            # nnpsav is the number of average power spectra in the file.
            nnpsav = t.size // avnn
            for n in range(nnpsav):
                nst0 = n * avnn
                nst1 = (n + 1) * avnn
                if (nst1 >= t.size):
                    nst1 = t.size
                taps = np.append(taps, (t[nst0] + t[nst1-1]) / 2.0)
                datet = np.append(datet, datet0 + \
                    np.timedelta64(np.int64(taps[-1] * 1.e6), 'us'))
                psbsum = np.append(psbsum, np.sum(sg[nf0:nf1+1, \
                    nst0:nst1].flatten()) / (nst1 - nst0))
    if (len(sgfiles) > 1):
        # Process the last file.
        info = get_sginfo(sgfiles[-1])
        if (type(info['datestring']) == np.bytes_) or \
            (type(info['datestring']) == bytes):
            info['datestring'] = info['datestring'].decode('utf-8')
        if (info['datestring'][-1] == 'Z'):
            if (not utc):
                warnings.warn('Inconsitent UTC time encoding') 
            # datestring is in UTC. Avoid depreciation warning in np.datetime64.
            datet0 = np.datetime64(info['datestring'][:-1])
        else:
            if (utc):
                warnings.warn('Inconsitent UTC time encoding') 
            # datestring is in local time.
            datet0 = np.datetime64(info['datestring'])
        nskip = info['nfft'] - info['noverlap']
        dt = nskip / info['fs']
        f, t, sg = sgread(sgfiles[-1], t0=0, t1=t1)
        t0f = t[0] - dt
        # nf0 is the index of the lowest frequency band to be included in the sum.
        nf0 = np.int64(np.floor(f0 * info['nfft'] / info['fs']))
        # nf1 is the index of the highest frequency band to be included in the sum.
        nf1 = np.int64(np.floor(f1 * info['nfft'] / info['fs']))
        # Determine the number of average power spectra in for the file.
        if (avtime == -1):
            taps = np.append(taps, (t[0] + t[-1]) / 2.0)
            datet = np.append(datet, datet0 + \
                np.timedelta64(np.int64(taps[-1] * 1.e6), 'us'))
            psbsum = np.append(psbsum, \
                np.sum(sg[nf0:nf1+1].flatten()) / t.size)
        else:
            # Find the number of spectrograms for which the last spectrogram is based on
            # samples entirely within the avtime window.
            avnn = np.int64(np.floor((avtime * info['fs'] - info['nfft']) / 
                (info['nfft'] - info['noverlap']) + 1))
            # nnpsav is the number of average power spectra in the file.
            nnpsav = t.size // avnn
            for n in range(nnpsav):
                nst0 = n * avnn
                nst1 = (n + 1) * avnn
                if (nst1 >= t.size):
                    nst1 = t.size
                taps = np.append(taps, (t[nst0] + t[nst1-1]) / 2.0)
                datet = np.append(datet, datet0 + \
                    np.timedelta64(np.int64(taps[-1] * 1.e6), 'us'))
                psbsum = np.append(psbsum, np.sum(sg[nf0:nf1+1, \
                    nst0:nst1].flatten()) / (nst1 - nst0))
    if savefile:
        if (foutname==-1):
            from os.path import basename
            foutname = basename(sgfiles[0]).rsplit('_sg', 1)[0] + \
                'psb_' + str(int(f0)) + '_' + str(int(f1)) + '.csv'
        psbsum_save(foutname, psbsum, f0, f1, taps, datet, utc, info)
    return(taps, datet, utc, psbsum)

def psbsum_save(fname, psbsum, f0, f1, t, datet, utc=True, sginfo=False, dB=True):
    import numpy as np
    # Make an array of datestrings with ms precision.
    datestr = np.array(np.array(datet, dtype='datetime64[ms]'), dtype=str)
    # Make an output array with rows for each row in the CSV file.
    output = np.zeros((psbsum.size,), dtype=[('datestring', '|S42'), \
        ('t', 'float64'), ('psbsum', 'float64')])
    output['datestring'] = datestr
    output['t'] = t
    if dB:
        output['psbsum'] = 10. * np.log10(psbsum)
    else:
        output['psbsum'] = psbsum   
    header = 'f0: ' + str(f0) + ', f1: ' + str(f1) + '\n' 
    if utc:
        header = header + 'Date String (UTC), Time (s)'
    else:
        header = header + 'Date String (local time), Time (s)'        
    if (dB):
        header = header + ', PSB sum (dB)'
    else:
        header = header + ', PSB sum'
    np.savetxt(fname, output, fmt='%s, %.18g, %.18g', delimiter=',', header=header)

def psbsum_plot_files(sgfiles, f0=0, f1=-1, t0=0, t1=-1, avtime=-1, dB=True,
    dtformat=None, savefile=False, foutname=None):
    """Use psbsum to produce a power spectral band sum and plot it.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from datetime import datetime
    from os.path import basename, dirname
    import warnings
    taps, datet, utc, psb = psbsum(sgfiles, f0, f1, t0, t1, avtime)
    # Need to convert from np.datetime64 to datetime.datetime to plot
    # with matplotlib.
    dtdatet = np.array([])
    for adatet in datet:
        dtdatet = np.append(dtdatet, _dt_convert(adatet))
    fig, ax = plt.subplots()
    if (dB):
        psbplot = plt.plot(dtdatet, 10*np.log10(psb), '-k')
        ax.set_ylabel('Power Spectral Density (dB)')
    else:
        psbplot = plt.plot(dtdatet, psb, '-k')
        ax.set_ylabel('Power Spectral Density')
    ax.xaxis_date()
    if (dtformat):
        date_format = mdates.DateFormatter(dtformat)
    else:
        if ((datet[-1] - datet[0]).astype('m8[D]') >
            np.timedelta64(0,'D')):
            # There is more than one day of data. Include date.
            date_format = mdates.DateFormatter('%Y-%m-%d %H:%M')
        else:
            date_format = mdates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(date_format)
    if (utc):
        ax.set_xlabel('UTC Time')
    else:
        ax.set_xlabel('Time')
    if (f1 == -1):
        ax.set_title('Frequency: > {0:.6g} Hz'.format(f0))
    else:
        ax.set_title('Frequency: {0:.6g} Hz - {1:.6g} Hz'.format(f0,
            f1)) 
    fig.autofmt_xdate()
    if (savefile):
        if (not foutname):
            foutname = basename(sgfiles[0]).rsplit('_sg', 1)[0] + \
                'psb_' + str(int(f0)) + '_' + str(int(f1)) + '.pdf'
        fig.savefig(foutname)
        plt.close(fig)
    else:
        fig.show()
        return(fig, ax, psbplot)

def multi_psbsum(sgfiles, fbands=[[0, -1]], t0=0, t1=-1, avtime=-1,
    total=True, dB=False, savefile=True, foutname=None):
    """Use psbsum to compute multiple power spectral band sums from the
    same input files.
    """
    import numpy as np
    from datetime import datetime
    from os.path import basename, dirname
    if isinstance(sgfiles, str):
        # There is only one file not in a list. Put it in a list.
        sgfiles = [sgfiles]
    elif not (isinstance(sgfiles, list) or isinstance(sgfiles,
        tuple)):
        raise ValueError('Parameter not list, tuple, or string',
            sgfiles)
    fbnp = np.array(fbands)
    if (fbnp.size == 2):
        # There is only one band entered as a 1-D list.
        fbnp = np.array(fbnp).reshape((1,2))
    if (fbnp.shape[1] != 2):
        raise ValueError('Parameter not a list of 2-value lists',
            fbands)
    if (total):
        fbnp = np.append(fbnp, np.array([[0, -1]]), axis=0)
    nnfb = fbnp.shape[0]
    taps, datet, utc, psb = psbsum(sgfiles, fbnp[0,0], fbnp[0,1], t0,
        t1, avtime)
    multipsb = np.array([psb])
    for nfb in range(1, nnfb):
        _, _, _, psb = psbsum(sgfiles, fbnp[nfb,0], fbnp[nfb,1], t0,
            t1, avtime)
        multipsb = np.append(multipsb, np.array([psb]), axis=0)
    if (savefile):
        if (not foutname):
            if (dirname(sgfiles[0]) == ''):
                foutname = basename(sgfiles[0]).rsplit('_sg', 1)[0] + \
                    'psbsums.csv'
            else:
                foutname = dirname(sgfiles[0]) + '/' + \
                    basename(sgfiles[0]).rsplit('_sg', 1)[0] + \
                    'psbsums.csv'
        # Make an array of datestrings with ms precision.
        datestr = np.array(np.array(datet, dtype='datetime64[ms]'),
            dtype=str)
        # Make an output array with rows for each row in the CSV file.
        dtp = [('datestring', '|S42'), ('t', 'float64')]
        for nfb in range(nnfb):
            dtp = dtp + [('psbsum[{}]'.format(nfb),'<f8')]
        dtp=np.dtype(dtp)
        output = np.zeros((psb.size,), dtype=dtp)
        output['datestring'] = datestr
        output['t'] = taps
        fmt = '%s, %.18g'
        if (dB):
            multipsb = 10. * np.log10(multipsb)
            psbunit= ' (dB)'
        else:
            psbunit = ''
        if utc:
            header = 'Date String (UTC), Time (s)'
        else:
            header = 'Date String (local time), Time (s)'        
        for nfb in range(nnfb):
            if (fbnp[nfb,1] == -1):
                if (fbnp[nfb,0] == 0):
                    header = header + ', Full Spectrum PSB Sum' + \
                        psbunit
                else:
                    header + (', PSB Sum f > {0:.6g} Hz' + \
                psbunit).format(fbnp[nfb,0])
            else:
                header = header + \
                    (', PSB Sum {0:.6g} Hz - {1:.6g} Hz' + \
                    psbunit).format(fbnp[nfb,0], fbnp[nfb,1])
            fmt = fmt + ', %.18g'
            psbkey = 'psbsum[{}]'.format(nfb)
            output[psbkey] = multipsb[nfb]
        np.savetxt(foutname, output, fmt=fmt, delimiter=',',
            header=header)
    else:
        return(taps, datet, utc, multipsb)
        


def _dt_convert(dt64):
    """Convert np.datetime64 ns (dtype('<M8[ns]') object to
    datetime.datettime object. 
    """
    from datetime import datetime
    ns = 1.e-9
    return(datetime.utcfromtimestamp(dt64.astype(int) * ns))
