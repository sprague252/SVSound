"""Module containing get_info program for a Broadcast Wave file written
by an AudioMoth device.
"""
from __future__ import division

def get_info(file, info={}):
    """info = get_info(file, info={})
    Read the information in a WAV file written by an AudioMoth recorder
    and return the contents. The standard information in the fmt chunk
    is included in the info dictionary along with other information
    encoded in INFO chunk.
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
        "Nsamples" - integer number of samples in the file (in each channel)
        "wavetype" - string with "AudioMoth" as the wave file type read.
        
    Additional info dictionary keys and values returned
    
        "ICMT" - string with the contents of the ICMT subchunk.
        "IART" - string with the contents of the IART subchunk.
        "datestring" - string with the date and time of the beginning
            of the recording in ISO 8601 format.
        "voltage" - string with the battery voltage at the beginning
            of the recording.
        "gain" - string with the AudioMoth gain setting for
            recording.
        "serial number" - string with the serial number of the
            AudioMoth recording device.      
    """
    
    import struct
    import numpy as np
    import datetime as dt
    import re
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
    if data != b'fmt ':
        raise ValueError('Chunk ID not fmt : ', data)
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
    if data != b'LIST':
        raise ValueError('Chunk ID not LIST: ', data)
    data = file.read(4)
    info['LISTsize'] = struct.unpack('<I',data)[0]
    data = file.read(4)
    if data != b'INFO':
        raise ValueError('Chunk Format not INFO: ', data)
    # The LIST chunk consists of INFO containing ICMT and IART
    # subchunks.
    data = file.read(4)
    while data != b'data':
        if data == b'ICMT':
            # This is the information subchunk. Store in a string.
            size = struct.unpack('<I', file.read(4))[0]
            info['ICMT'] = (
                file.read(size).decode('ascii').rstrip('\x00'))
            # Extract the recording time.
            pattern = (r'(\d\d):(\d\d):(\d\d) ' +
                r'(\d\d)/(\d\d)/(\d\d\d\d) ' +
                r'\(UTC([+-])(\d\d?)\)')
            tmatch = re.search(pattern, info['ICMT'])
            if len(tmatch[8]) == 1:
                tz = tmatch[7] + '0' + tmatch[8] + '00'
            else:
                tz = tmatch[7] + tmatch[8] + '00'
            tstring = (tmatch[6] + tmatch[5] + tmatch[4] + 'T' +
                tmatch[1] + tmatch[2] + tmatch[3] + tz)
            info['datestring'] = dt.datetime.strptime(tstring,
                '%Y%m%dT%H%M%S%z').isoformat()
            # Extract the battery voltage.
            pattern = r'\d\.\dV'
            info['voltage'] = re.search(pattern, info['ICMT'])[0]
            # Extract the temperature.
            pattern = r'-?\d\d?\.\dC'
            info['temperature'] = re.search(pattern, info['ICMT'])[0]
            # Extract the gain setting.
            pattern = r'(\w+-?\w*) gain'
            info['gain'] = re.search(pattern, info['ICMT'])[1]
        elif data == b'IART':
            # This is the artist subchunk. Store in a string.
            size = struct.unpack('<I', file.read(4))[0]
            info['IART'] = (
                file.read(size).decode('ascii').rstrip('\x00'))
            # Extract the serial number
            pattern = r'AudioMoth ([0-9A-Z]+)'
            info['serial number'] = re.search(pattern, info['IART'])[1]
        data = file.read(4)
    # We are now at the data chunk.
    # Nsamples is the number samples in the file (should be total time * fs).
    info['Nsamples'] = struct.unpack('<I', file.read(4))[0] // \
        (info['block_align'])
    # data0 is the position (in bytes) of the first data sample.
    info['data0'] = file.tell()
    info['wavetype'] = 'AudioMoth'
    return(info)
    