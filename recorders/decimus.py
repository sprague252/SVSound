"""Package containing get_info program for a Broadcast Wave file written
by a Decimus device.
"""
from __future__ import division

def get_info(file, info={}):
    """Read the information in a Decimus WAV file, and return the contents.
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
    while data != b'fmt ':
        chunk_name = 'Chunk ' + data
        chunk_size = struct.unpack('<I', file.read(4))[0]
        data = file.read(chunk_size)
        info[chunk_name] = data
        data = file.read(4)
    fmt_size = struct.unpack('<I', file.read(4))[0]
    info['compress'] = struct.unpack('<H', file.read(2))[0]
    if info['compress'] != 1:
        raise ValueError('Only linear PCM files supported. Compress: ', info['compress'])
    info['chan'] = struct.unpack('<H', file.read(2))[0]
    info['fs'] = struct.unpack('<I', file.read(4))[0]
    info['byte_per_s'] = struct.unpack('<I', file.read(4))[0]
    info['block_align'] = struct.unpack('<H', file.read(2))[0]
    info['bits'] = struct.unpack('<H', file.read(2))[0]
    if fmt_size > 16:
        # Some devices add extra information in the fmt chunk. Read it and keep going.
        size = struct.unpack('<H', file.read(2))[0]
        info['extra fmt'] = file.read(size)
    data = file.read(4)
    if data == b'PAD ':
        # Some devices insert a padding chunk. Skip it.
        pad_size = struct.unpack('<I', file.read(4))[0]
        file.seek(pad_size,1)
        data = file.read(4)
    if data != b'data':
        raise ValueError('Chunk ID not data: ', data)
    # We are at the data chunk.
    # Nsamples is the number samples in the file (should be total time * fs).
    info['Nsamples'] = struct.unpack('<I', file.read(4))[0] // \
        (info['block_align'])
    # data0 is the position (in bytes) of the first data sample.
    info['data0'] = file.tell()
    info['wavetype'] = 'decimus'
    return(info)