"""Read and write data to/from Broadcast Wave files.  The following
propriatory boadcast wave file formats are currently supported:

  generic - generic Windows WAVE file format containing the basic
       infromation in the WAVE format chunk. No additional metadata is read.
       
  decimus - Windows WAVE file written by the Decimus(R) passive
      acoustics monitoring system and other devices that use the SA
      Instrumentation DAQ card. Metadata is extracted from the filename,
      which includes a timestamp, into the info dictionary.

  icListen - WAVE files written by icListen(R) recording devices.
      Metadata in the INFO chunk is read into the info dictionary.

  zoom - WAVE files written by ZOOM(R) recording devices. Metadata
      from the bext chunk and the iXML chunk is read into the info
      dictionary.
"""

from __future__ import division, print_function

def read(filename, t0=0, t1=-1, wavetype=None, chunk_b=3072):
    """info, wave = read(filename, t0, t1, wavetype, chunk_b)
    
    Read a WAV file and return the file information and waveform data.
    Input parameters
    
      filename - string with the name of the input WAV file
      t0 - start time in seconds for returned data (default: 0)
      t1 - end time in seconds for returned data. Value of -1 represents
        the end of the file. (default: -1)
      wavetype - string representing the type of WAV file (default:
        None). Currnetly supported types are 'generic', 'decimus',
        'icListen', and 'zoom'. If the value is None, the wavetype
        is determined using identify.
      chunk_b - number of bytes for each data chunk read from the
        file (default: 3072)
     
    Output
    
      info - dictionary with file information
      wave - Numpy array with waveform data values
    """
    import struct
    import numpy as np
    file = open(filename, 'rb')
    info = getInfo(file, wavetype)
    wave = wave_chunk(file, info, t0, t1, chunk_b=chunk_b)
    file.close()
    return(info, wave)

def wavread(filename, t0=0, t1=-1, wavetype=None, chunk_b=3072):
    """Depricated: This is a frontend to wavefile.read maintained for
    legacy purposes. See wavefile.read for complete usage information.
    """
    import warnings
    deprication_msg = "wavefile.wavread() is depricated. Use wavefile.read()."
    warnings.warn(deprication_msg, FutureWarning)
    info, wave = read(filename, t0, t1, wavetype, chunk_b)
    return(info, wave)


def identify(file):
    """wavetype = identify(file)
    
    Identify the type of WAV file and return its type. Files that are
    unable to be identified are classified as generic. The wave type
    identification allows the extraction of proprietary metadata stored
    in the file and filename.
    Input parameter
    
      file - filehandle for the WAV file to be identified
            
    Output

      wavetype - string with the name of the wave file type."""
    import warnings
    from os.path import basename
    import struct
    # Test each known file type for type-specific properties.
    # First make sure the file is a WAVE file.
    file.seek(0)
    data = file.read(4)
    if data != b'RIFF':
        if data[0:3] == b'ID3':
            # There is an ID3 tag at the beginning if the file.
            # Only ID3v2 tags are supported. See <https://id3.org/id3v2.3.0>.
            file.seek(2, 1)
            data0 = file.read(1)
            data1 = file.read(1)
            data2 = file.read(1)
            data3 = file.read(1)
            tag_size = (
                int.from_bytes(data3, 'big') + int.from_bytes(data2,
                'big') * 2**7 + int.from_bytes(data1, 'big') * 2**14 +
                int.from_bytes(data0, 'big') * 2**21
           )
            # Skip forward tag_size bytes.
            file.seek(tag_size, 1)
            data = file.read(4)
            if data != b'RIFF':
                raise ValueError('Chunk ID not RIFF: ', data)
        else:
            raise ValueError('Chunk ID not RIFF: ', data)
    file.seek(4, 1)
    data = file.read(4)
    if data != b'WAVE':
        raise ValueError('Chunk ID not WAVE: ', data)
    data = file.read(4)
    # Check to see if the file was written by an icListen device.
    if data == b'LIST':
        # Looks like an icListen file.  Check to see if it has an INFO chunk.
        file.seek(4, 1)
        data = file.read(4)
        if data == b'INFO':
            wavetype = 'icListen'
        else:
            message = 'Not icListen...Chunk Format not INFO...Using generic format'
            warnings.warn(message)
            wavetype = 'generic'
    # Check to see if the this is a bext chunk (used by ZOOM devices).
    elif data == b'bext':
        file.seek(260, 1)
        data = file.read(32)
        tstr = data.split(b'\0',1)[0].split(b' ')[0]
        if tstr == b'ZOOM':
            wavetype = 'zoom'
        else:
            message = 'bext originator string not ZOOM...Using generic format'
            warnings.warn(message)
            wavetype = 'generic'
    # Check to see if it is a Decimus file. Decimus files have names
    # beginning with 'PAM'
    elif basename(file.name)[0:3] == 'PAM':
        wavetype = 'decimus'
    else:
        # Use generic WAVE file format.
        wavetype = 'generic'
    return(wavetype)

def getInfo(file, wavetype=None):
    """Read the file information from the WAV file with handle file and
    return the info dictionary.
    """
    if wavetype == None:
        wavetype = identify(file)
    info = {'wavetype' : wavetype}
    if wavetype == 'decimus':
        from .recorders.decimus import get_info
    if wavetype == 'icListen':
        from .recorders.icListen import get_info
    if wavetype == 'zoom':
        from .recorders.zoom import get_info
    if wavetype == 'generic':
        from .recorders.generic import get_info
    info = get_info(file, info)
    return(info)
    
def wave_chunk(file, info, t0=0, t1=-1, chunk_b=768, verbose=False):
    """Read a WAVE file in chunks (not all at once) and return all the
    data. This is a back-end to the read function.
    """
    import numpy as np
    import struct
    import logging
    ch = logging.StreamHandler()
    formatter = logging.Formatter('%(asctime)s %(levelname)8s %(name)s | %(message)s')
    ch.setFormatter(formatter)
    logger = logging.getLogger('wave_chunk')
    logger.addHandler(ch)
    if verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)
    t0_pos = info['data0'] + np.uint32(np.floor(t0 * info['byte_per_s']))
    if t1==-1:
        t1_pos = info['data0'] + info['Nsamples'] * info['block_align']
    else:
        t1_pos = info['data0'] + np.uint32(np.floor(t1 * info['byte_per_s'])) + 1
    if (t1_pos > info['filesize']):
        try:
            raise ValueError('Past EOF')
        except ValueError:
            print("t1 is past the end of file")
            print("t1_pos: ", t1_pos)
            print("info['filesize']: ", info['filesize'])
            raise
    # nn is the number of bytes to read.
    nn = t1_pos - t0_pos
    logger.info('wave_chunk reading ' + str(nn) + 'bytes')
    #nnchunk is the number of full chunks to read.
    nnchunk = nn // chunk_b
    logger.info(str(nnchunk) + ' chunks')
    # ncrumb is the leftovers that did not fit in the last chunk.
    ncrumb = nn % chunk_b
    fmt0 = '<' + str(chunk_b * 8 // info['bits'])
    # Determine how to decode data.
    if info['bits'] == 8:
        wave = np.array([], dtype = np.int8)
        decode = lambda bindat: struct.unpack(fmt0 + 'b', bindat)
    elif info['bits'] == 16:
        wave = np.array([], dtype = np.int16)
        decode = lambda bindat: struct.unpack(fmt0 + 'h', bindat)
    elif info['bits'] == 24:
        wave = np.array([], dtype = np.int32)
        decode = _decode24
    elif info['bits'] == 32:
        if info['compress'] == 1:
            # 32 bit signed integer
            wave = np.array([], dtype = np.int32)
            decode = lambda bindat: struct.unpack(fmt0 + 'i', bindat)
        elif info['compress'] == 3:
            # 32 bit float
            wave = np.array([], dtype = np.float32)
            decode = lambda bindat: struct.unpack(fmt0 + 'f', bindat)
        else:
            raise ValueError('Not a recognized format. Bits: ', 
                info['bits'], ' Compress: ', info['compress'])
    elif info['bits'] == 64:
        if info['compress'] == 3:
            # 64 bit float
            wave = np.array([], dtype = np.float64)
            decode = lambda bindat: struct.unpack(fmt0 + 'd', bindat)
        else:
            raise ValueError('Not a recognized format. Bits: ', 
                info['bits'], ' Compress: ', info['compress'])
    else:
        raise ValueError('Not a recognized format. Bits: ', 
            info['bits'], ' Compress: ', info['compress'])
    # Now start reading the data.
    file.seek(t0_pos)
    for nchunk in range(nnchunk):
        logger.info('  Chunk ' + str(nchunk) + ' of ' + str(nnchunk))
        chunk = file.read(chunk_b)
        wave = np.append(wave, decode(chunk))
    if ncrumb:
        # Check to see that ncrumb is a whole number of samples
        if (np.mod(ncrumb, info['bits']//8) >= 0):
            # ncrumb has a partial sample at the end likely due to truncation error.
            logger.info('Partial sample in ncrumb:')
            logger.info('  {:d} bytes out of {:d} bytes'.format(np.mod(ncrumb, \
                info['bits']//8), info['bits']//8))
            if (t1 == -1 or t1 >= (info['Nsamples'] - 1) / info['fs']):
                # The last sample in the crumb is the last sample in the file.
                # Reduce ncrumb to use the last complete sample.
                ncrumb = np.int32(np.floor(ncrumb / info['bits']/8)) * \
                    info['bits']//8
                logger.info('Truncating ... ncrumb now {:d}'.format(ncrumb))
            else:
                # There is at least one more sample available in the
                # file. Use it to complete ncrumb.
                ncrumb = np.int32(np.ceil(ncrumb / info['bits']/8)) * \
                    info['bits']//8
                logger.info('Extending ... ncrumb now {:d}'.format(ncrumb))
        logger.info('  Cleaning up ' + str(ncrumb) + ' bytes')
        fmt0 = '<' + str(ncrumb * 8 // info['bits'])
        crumb = file.read(ncrumb)
        wave = np.append(wave, decode(crumb))
    logger.info('wave_chunk done')
    if info['chan'] > 1:
        # Multichannel file. Consecutive samples are in different
        # channels. Reshape with channels in different rows.
        logger.info('Sorting channels')
        wave = np.reshape(wave, (info['chan'], -1), order='F')
    logger.removeHandler(ch)
    return(wave)

def _decode24(string):
    """Decode 24-bit binary integer data, often used in WAV files."""
    from struct import unpack
    fmt = '<'+ str(len(string)//3) + 'i'
    # Convert the string to a 32 bit (4-byte) string by adding \0 as the
    # most significant byte for positive numbers and \xff as the most
    # significant byte for negative numbers.
    str32=b''
    for n in range(len(string))[:-1:3]:
        str32 += string[n:n+3] + (b'\0' if string[n+2] < 128 else b'\xff')
    # Now use the 32-bit unpack function ('i' format) to get the number.
    return(unpack(fmt, str32))
