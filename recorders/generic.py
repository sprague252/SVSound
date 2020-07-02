"""Module containing get_info program for a generic Broadcast Wave file.
"""
from __future__ import division

def get_info(file, info={}):
    """info = get_info(file, info={})
    Read the information in a generic WAV file, and return the
    contents. Only the standard information in the fmt chunk is included
    in the info dictionary.
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
        "wavetype" - string with "generic" as the wave file type read.
    """
    import struct
    import numpy as np
    # tag_size is for ID3 tags
    tag_size = None
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
    data = file.read(4)
    info['filesize'] = struct.unpack('<I',data)[0] + 8
    if tag_size != None:
        info['filesize'] = info['filesize'] + tag_size + 10
    data = file.read(4)
    if data != b'WAVE':
        raise ValueError('Chunk ID not WAVE: ', data)
    data = file.read(4)
    while data != b'fmt ':
        chunk_name = b'Chunk ' + data
        chunk_size = struct.unpack('<I', file.read(4))[0]
        data = file.read(chunk_size)
        info[chunk_name] = data
        data = file.read(4)
    fmt_size = struct.unpack('<I', file.read(4))[0]
    info['compress'] = struct.unpack('<H', file.read(2))[0]
    if (info['compress'] != 1 and info['compress'] != 3):
        raise ValueError('Only integer and floating point linear PCM ' +
        'files supported. Compress: ', info['compress'])
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
    if (data == b'PAD ' or data == b'FLLR'):
        # Some devices insert a padding chunk or a filler chunk. Skip it.
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
    info['wavetype'] = 'generic'
    return(info)
