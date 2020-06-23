"""Package containing get_info program for a Broadcast Wave file written
by an icListen device.
"""
from __future__ import division

def get_info(file, info={}):
    """Read the information in a Broadcast Wave Format (BWF) file
    written by the icListen recorder and return the contents.
    """
    import struct
    import numpy as np
    file.seek(0)
    data = file.read(4)
    if data != b'RIFF':
        raise ValueError('Chunk ID not RIFF: ', data)
    data = file.read(4)
    info['filesize'] = struct.unpack('<I',data)[0] + 8
    data = file.read(4)
    if data != b'WAVE':
        raise ValueError('Chunk ID not WAVE: ', data)
    data = file.read(4)
    if data != b'LIST':
        raise ValueError('Chunk ID not LIST: ', data)
    data = file.read(4)
    info['LISTsize'] = info['filesize'] = struct.unpack('<I',data)[0]
    data = file.read(4)
    if data != b'INFO':
        raise ValueError('Chunk Format not INFO: ', data)
    # Now start reading info chunks and storing them in info{}.
    key = file.read(4)
    while key != b'fmt ':
        size = struct.unpack('<I', file.read(4))[0]
        info[key] = file.read(size)
        key = file.read(4)
    # We are now in the fmt sub chunk. Set some INFO parameters before 
    # continuing.
    icmt = info['ICMT'].split()
    vmax = np.float64(icmt[0])
    # The following values are encoded as integers but can be forced into 
    # float64 (double precision) with no loss of precision.
    sens = np.float64(icmt[3])
    wavemax = np.float64(icmt[13])
    info['cal'] = vmax / wavemax * 10**(-sens / 20.0)
    # Continue with fmt sub chunk.
    fmtsize = struct.unpack('<I', file.read(4))[0]
    info['compress'] = struct.unpack('<H', file.read(2))[0]
    info['chan'] = struct.unpack('<H', file.read(2))[0]
    info['fs'] = struct.unpack('<I', file.read(4))[0]
    info['byte_per_s'] = struct.unpack('<I', file.read(4))[0]
    info['block_align'] = struct.unpack('<H', file.read(2))[0]
    info['bits'] = struct.unpack('<H', file.read(2))[0]
    if fmtsize > 16:
        size = struct.unpack('<H', file.read(2))[0]
        info['extra fmt'] = file.read(size)
    data = file.read(4)
    if data != b'data':
        raise ValueError('Chunk ID not data: ', data)
    # We are at the data chunk.
    # Nsamples is the number samples in the file (should be total time * fs).
    info['Nsamples'] = struct.unpack('<I', file.read(4))[0] // \
        (info['block_align'])
    # data0 is the position (in bytes) of the first data sample.
    info['data0'] = file.tell()
    info['wavetype'] = 'icListen'
    return(info)
