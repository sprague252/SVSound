"""Package containing get_info and wavread programs for a Broadcast Wave
file written by an ZOOM F8 recorder.
"""
from __future__ import division

def get_info(file, info={}):
    """Read the information in a Broadcast Wave Format (BWF) file
    written by the ZOOM F8 recorder and return the contents.
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
    if data != b'bext':
        raise ValueError('Chunk ID not bext: ', data)
    bext_size = struct.unpack('<I', file.read(4))[0]
    data = file.read(256)  
    info['desc'] = data.split(b'\0',1)[0]
    data = file.read(32)
    info['Originator'] = data.split(b'\0',1)[0]
    data = file.read(32)
    info['OriginatorReference'] = data.split(b'\0',1)[0]
    info['OriginationDate'] = file.read(10)
    info['OriginationTime'] = file.read(8)
    info['TimeReferenceLow'] = struct.unpack('<I', file.read(4))[0]
    info['TimeReferenceHigh'] = struct.unpack('<I', file.read(4))[0]
    info['Version'] = struct.unpack('<H',file.read(2))[0]
    info['UMID'] = file.read(64)
    info['LoudnessValue'] = struct.unpack('<H',file.read(2))[0]
    info['LoudnessRange'] = struct.unpack('<H',file.read(2))[0]
    info['MaxTruePeakLevel'] = struct.unpack('<H',file.read(2))[0]
    info['MaxMomentaryLoudness'] = struct.unpack('<H',file.read(2))[0]
    info['MaxShortTermLoudness'] = struct.unpack('<H',file.read(2))[0]
    file.seek(180,1)
    data = file.read(bext_size - 602)
    info['CodingHistory'] = data.split(b'\0',1)[0]
    data = file.read(4)
    if data != b'iXML':
        raise ValueError('Chunk ID not iXML: ', data)
    iXML_size = struct.unpack('<I', file.read(4))[0]
    data = file.read(iXML_size)
    info['iXML'] = data.split(b'\0',1)[0]
    data = file.read(4)
    if data != b'fmt ':
        raise ValueError('Chunk ID not fmt : ', data)
    fmt_size = struct.unpack('<I', file.read(4))[0]
    info['compress'] = struct.unpack('<H', file.read(2))[0]
    info['chan'] = struct.unpack('<H', file.read(2))[0]
    info['fs'] = struct.unpack('<I', file.read(4))[0]
    info['byte_per_s'] = struct.unpack('<I', file.read(4))[0]
    info['block_align'] = struct.unpack('<H', file.read(2))[0]
    info['bits'] = struct.unpack('<H', file.read(2))[0]
    if fmt_size > 16:
        size = struct.unpack('<H', file.read(2))[0]
        info['extra fmt'] = file.read(size)
    data = file.read(4)
    if data == b'PAD ':
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
    info['wavetype'] = 'zoom'
    return(info)