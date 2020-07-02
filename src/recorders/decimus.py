"""Module containing get_info program for a Broadcast Wave file written
by a Decimus device."""
from __future__ import division

def get_info(file, info={}):
    """info = get_info(file, info={})
    Read the information in a WAV file written by a Decimus device and
   return the contents. The standard information in the fmt chunk is
   included in the info dictionary.
    Input Parameters
    
        file - filehandle of an open WAV file
        info - (optional) dictionary that may contain file
            information from other sources. Defaults to an empty
            dictionary.
    
    Output
    
        info - dictionary with information read from the file. If an
            info dictionary was supplied as an input parameter,
            entires that were not changed are also included.
    
    info dictionary keys and values returned
    
        "bits" - integer with the number of bits in each sample.
        "block_align" - number of bytes sampled at the same time (all
            channels combined) in the data
        "byte_per_s" - integer number of bytes per second recorded
        "chan" - integer number of channels in the file
        "compress" - integer Wave file compression index. Only 1
            (uncompressed integer data) and 3 (uncompressed floating
            point data) are currently supported.
        "data0" - integer byte address of the first sample in the file
        "filesize" - integer size of the file in bytes
        "fs" - integer sample rate in samples/second
        "Nsamples" - integer number of samples in the file (in each channel)
        "wavetype" - string with "decimus" as the wave file type read.
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