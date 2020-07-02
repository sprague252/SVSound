"""Module containing get_info and wavread programs for a Broadcast Wave
file written by an ZOOM F8 recorder.
"""
from __future__ import division

def get_info(file, info={}):
    """Read the information in a Broadcast Wave Format (BWF) file
    written by the ZOOM F8 recorder and return the contents.
    """
    """info = get_info(file, info={}) Read the information in a WAV file
    written by a Zoom device (this function was written specifically for
    recordings made by the Zoom F8 recorder but it may work with other
    Zoon devices as well) and return the contents. The standard
    information in the fmt chunk is included in the info dictionary as
    well as information in the bext chunk and the iXML chunk.
    Input Parameters
    
        file - filehandle of an open WAV file
        info - (optional) dictionary that may contain file
            information from other sources. Defaults to an empty
            dictionary.
    
    Output
    
        info - dictionary with information read from the file. If an
            info dictionary was supplied as an input parameter,
            entires that were not changed are also included.
    
    Standard info dictionary keys and values returned
    
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
        "Nsamples" - integer number of samples in the file (in each
            channel)
        "wavetype" - string with "decimus" as the wave file type read.
    
    bext Chunk info dictionary keys and values returned. (See Zoom
    documentation for details.)
    
        "CodingHistory" - coding history string
        "desc" - recording description string
        "LoudnessRange" - int16 recording loudness range value
        "LoudnessValue" - int16 recording loudness value
        "MaxMomentaryLoudness" - int16 recording maximum momentary
            loudness value
        "MaxShortTermLoudness" - int16 recording maximum short term
            loudness value
        "MaxTruePeakLevel" - int16 recording maximum maximum true
            peak level
        "OriginationDate" - recording origination date string
        "OriginationTime" - recording origination time string
        "Originator" - recording originator string
        "OriginatorReference" - recording originator reference string
        "TimeReferenceHigh" - int32 time of high sample in recording
        "TimeReferenceLow" - int32 time of low sample in recording
        "UMID" - UMID string
        "Version" - int16 version number
        
    The contents in the entire iXML block are stored in info["iXML"] as
    a string."""

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