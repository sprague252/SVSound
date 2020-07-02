# SVSound

This is Python package for reading Broadcast Wave Files in various formats along with metadata written by several recording devices. The content in this module was forked from Mark Sprague's collection of sound recording and analysis modules and intended for use by his students. This package is available to everyone under the GPL 3.0 license. 

## wavefile Module

The module `wavefile` contains programs for reading Broadcast Wave (.wav) files.
The following propriatory boadcast wave file formats are currently
supported:

* generic - generic Windows WAVE file format containing the basic
infromation in the WAVE format chunk. No additional metadata is read.

* decimus - Windows WAVE file written by the Decimus(R) passive
acoustics monitoring system and other devices that use the SA
Instrumentation DAQ card. Metadata is extracted from the filename, which
includes a timestamp, into the info dictionary.

* icListen - WAVE files written by icListen(R) recording devices.
Metadata in the INFO chunk is read into the info dictionary.

* zoom - WAVE files written by ZOOM(R) recording devices. Metadata from
the bext chunk and the iXML chunk is read into the info dictionary.

### Functions

#### read( )

`info, wave = read(filename, t0, t1, wavetype, chunk_b)`

Read a WAV file and return the file information and waveform data. This function includes support for single and multiple channel files encoded in linear PCM format with the following data formats (all little-endian):

  * 8 bit signed integer
  * 16 bit signed integer
  * 24 bit signed integer
  * 32 bit signed integer
  * 32 bit floating point
  * 64 bit floating point

Input parameters
    
>`filename` - string with the name of the input WAV file

>`t0` - start time in seconds for returned data (default: 0)

>`t1` - end time in seconds for returned data. Value of -1 represents
        the end of the file. (default: -1)
        
>`wavetype` - string representing the type of WAV file (default:
None). Currnetly supported types are 'generic', 'decimus',
'icListen', and 'zoom'. If the value is None, the wavetype
is determined using identify.

>`chunk_b` - number of bytes for each data chunk read from the
file (default: 3072)
     
Output
    
>`info` - dictionary with file information

>`wave` - Numpy array with waveform data values

#### identify( )

`wavetype = identify(file)`
    
Identify the type of WAV file and return its type. Files that are unable
to be identified are classified as generic. The wave type identification
allows the extraction of proprietary metadata stored in the file and
filename.

Input parameter
    
>`file` - filehandle for the WAV file to be identified
            
Output

>`wavetype` - string with the name of the wave file type.

#### wave_chunk( )

Read a WAVE file in chunks (not all at once) and return all the data.
This is a back-end to the read function and is not intended for high-level use.

## recorders Subpackage

The subpackage `recorders` contains modules with specific `get_info()` functions for each supported recorder type. Currently supported recorders are described in the wavefile Module introduction (above). Each `get_info()` function has the same input and output parameters and usage.

`info = get_info(file, info)`

Read the information in a generic WAV file, and return the contents.
Only the standard information in the fmt chunk is included in the info
dictionary.

Input Parameters
    
>`file` - filehandle of an open WAV file

>`info` - (optional) dictionary that may contain file
    information from other sources. Defaults to an empty
    dictionary.
    
Output
    
>`info` - dictionary with information read from the file. If an info
dictionary was supplied as an input parameter, entires that were not
changed are also included.
    
Standard `info` dictionary keys and values returned for all file types:
    
>`"bits"` - integer with the number of bits in each sample.

>`"block_align"` - number of bytes sampled at the same time (all
channels combined) in the data

>`"byte_per_s"` - integer number of bytes per second recorded

>`"chan"` - integer number of channels in the file

>`"compress"` - integer Wave file compression index. Only 1
(uncompressed integer data) and 3 (uncompressed floating point data) are
currently supported.

>`"data0"` - integer byte address of the first sample in the file

>`"filesize"` - integer size of the file in bytes

>`"fs"` - integer sample rate in samples/second

>`"Nsamples"` - integer number of samples in the file (in each channel)

>`"wavetype"` - string with the name file type read.

Other keys and values in the `info` dictionary are recorder-specific and
depend on the `wavetype` value.

### Recorder-Specific info keys and values

#### Decimus

Recordings identified as Decimus recordings have `info["wavetype"]` set
to "decimus". Otherwise, they contain only the standard info keys and
values.

#### Generic

Recordings classified as generic have `info["wavetype"]` set to
"generic" and contain only the standard info keys and values.

#### icListen

Recordings identified as icListen recordings have `info["wavetype"]` set
to "icListen". In addition, each key/value pair written to the INFO
chunk in the file are added to `info`. See the icListen documentation
for details on these parameters. 

The value `info["cal"]` contains a float64 calibration value for the data.
Multiply data samples by this value to obtain calibrated values in
micropascals.

#### Zoom

Recordings identified as Zoom recordings have `info["wavetype"]` set to
"zoom". The following information encoded in the bext chunk is added to
`info` as keys and values (See Zoom documentation for details.)
    
>`"CodingHistory"` - coding history string

>`"desc"` - recording description string

>`"LoudnessRange"` - int16 recording loudness range value

>`"LoudnessValue"` - int16 recording loudness value

>`"MaxMomentaryLoudness"` - int16 recording maximum momentary
loudness value

>`"MaxShortTermLoudness"` - int16 recording maximum short term
loudness value

>`"MaxTruePeakLevel"` - int16 recording maximum maximum true
peak level

>`"OriginationDate"` - recording origination date string

>`"OriginationTime"` - recording origination time string

>`"Originator"` - recording originator string

>`"OriginatorReference"` - recording originator reference string

>`"TimeReferenceHigh"` - int32 time of high sample in recording

>`"TimeReferenceLow"` - int32 time of low sample in recording

>`"UMID"` - UMID string

>`"Version"` - int16 version number

The contents in the entire iXML block are stored in `info["iXML"]` as a
string.


