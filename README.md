# SVSound

This is Python package for reading Broadcast Wave Files in various formats along with metadata written by several recording devices. The content in this module was forked from Mark Sprague's collection of sound recording and analysis modules and intended for use by his students. This package is available to everyone under the GPL 3.0 license. 

## Versions

* 0.0.1 - initial release
* 0.0.2 - Corrected a bug that caused the end of a wave file to be truncated.
* 0.0.3 - Corrected a bug causing incorrect partial file read; include ID3 contents (in binary form) in info.
* 0.0.4 - Now skips over nonstandard file chunks without producing errors.
* 0.0.5 - contains an error ... do not use
* 0.0.6 - Not many changes other than support for Python 3.9.
* 0.0.7 - Support for Python 3.10; added levels subpackage to calculate SPL values and SEL values from sampled sounds.
* 0.0.8--0.0.10 - All contain an error ... do not use
* 0.1.0 -  Decimus, icListen, and Zoom recorders fully supported. Correctly converts binary strings to ACII stings when interpreting nonstandard sub-chunks written by icListen and Zoom recorders. Determines file length correctly. Module `levels` imports module `wavefile` correctly. Supports Python <= 3.11.
* 0.2.0 - AudioMoth recorders fully supported with data extracted from the ICMT and IART sub-chunks added to the info dictionary. Supports Python <= 3.12.
* 0.2.1 - Bumped version to fix errors in pypi upload.
* 0.3.0 - Added sgwave to do some spectral analysis.

## Installation

### Conda/Mamba

Install using `conda`:

```
conda install -c sprague252 svsound
```

Install using `mamba`:

```
mamba install -c sprague252 svsound
```

### Pip

Install using `pip`:

```
pip install svsound
```


## wavefile Module

The module `wavefile` contains programs for reading Broadcast Wave
(.wav) files. The following propriatory boadcast wave file formats are
currently supported:

* generic - generic Windows WAVE file format containing the basic
infromation in the WAVE format chunk. No additional metadata is read.

* AudioMoth - Windows WAVE file written by an AudioMoth recording
device. Metadata extracted from the ICMT and IART chunks is added to the
info dictionary.

* decimus - Windows WAVE file written by the Decimus&reg; passive
acoustics monitoring system and other devices that use the SA
Instrumentation DAQ card. Metadata is extracted from the filename, which
includes a timestamp, into the info dictionary.

* icListen - WAVE files written by icListen&reg; recording devices.
Metadata in the INFO chunk is read into the info dictionary.

* zoom - WAVE files written by ZOOM&reg; recording devices. Metadata from
the bext chunk and the iXML chunk is read into the info dictionary.

### Functions

#### read( )

`info, wave = read(filename, t0, t1, wavetype, chunk_b, verbose)`

Read a WAV file and return the file information and waveform data. This
function includes support for single and multiple channel files encoded
in linear PCM format with the following data formats (all
little-endian):

  * 8 bit signed integer
  * 16 bit signed integer
  * 24 bit signed integer
  * 32 bit signed integer
  * 32 bit floating point
  * 64 bit floating point

Input parameters
    
<p style="margin-left: 3em; text-indent: -2em;">
<code>filename</code> - string with the name of the input WAV file 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>t0</code> - start time in seconds for returned data (default: 0) 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>t1</code> - end time in seconds for returned data. Value of -1
represents the end of the file. (default: -1) 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>wavetype</code> - string representing the type of WAV file
(default: None). Currnetly supported types are 'generic', 'AudioMoth',
'decimus', 'icListen', and 'zoom'. If the value is None, the wavetype is
determined using identify. 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>chunk_b</code> - number of bytes for each data chunk read from the
file (default: 3072) 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>verbose</code> - give verbose status updates (default: False)
</p>


Output
    
<p style="margin-left: 3em; text-indent: -2em;">
<code>info</code> - dictionary with file information and metadata (if available)
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>wave</code> - Numpy array with waveform data values. For a single
channel file, <code>wave</code> is a flat, 1-D array. For a multichannel
recording each channel is a row in <code>wave</code>, so
<code>wave[0]</code> is the first channel, <code>wave[1]</code> the
second channel, etc.
</p>

#### identify( )

`wavetype = identify(file)`
    
Identify the type of WAV file and return its type. Files that are unable
to be identified are classified as generic. The wave type identification
allows the extraction of proprietary metadata stored in the file and
filename.

Input parameter
    
<p style="margin-left: 3em; text-indent: -2em;">
<code>file</code> - filehandle for the WAV file to be identified 
</p>
            
Output

<p style="margin-left: 3em; text-indent: -2em;">
<code>wavetype</code> - string with the name of the wave file type. 
</p>

#### wave_chunk( )

Read a WAVE file in chunks (not all at once) and return all the data.
This is a back-end to the read function and is not intended for
high-level use.

## recorders Subpackage

The subpackage `recorders` contains modules with specific `get_info()`
functions for each supported recorder type. Currently supported
recorders are described in the wavefile Module introduction (above).
Each `get_info()` function has the same input and output parameters and
usage.

`info = get_info(file, info)`

Read the information in a generic WAV file, and return the contents.
Only the standard information in the fmt chunk is included in the info
dictionary.

Input Parameters
    
<p style="margin-left: 3em; text-indent: -2em;">
<code>file</code> - filehandle of an open WAV file 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>info</code> - (optional) dictionary that may contain file
information from other sources. Defaults to an empty dictionary. 
</p>
    
Output
    
<p style="margin-left: 3em; text-indent: -2em;">
<code>info</code> - dictionary with information read from the file. If
an info dictionary was supplied as an input parameter, entires that were
not changed are also included. 
</p>
    
Standard `info` dictionary keys and values returned for all file types:
    
<p style="margin-left: 3em; text-indent: -2em;">
<code>"bits"</code> - integer with the number of bits in each sample. 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"block_align"</code> - number of bytes sampled at the same time
(all channels combined) in the data 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"byte_per_s"</code> - integer number of bytes per second recorded 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"chan"</code> - integer number of channels in the file 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"compress"</code> - integer Wave file compression index. Only 1
(uncompressed integer data) and 3 (uncompressed floating point data) are
currently supported. 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"data0"</code> - integer byte address of the first sample in the
file 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"filesize"</code> - integer size of the file in bytes 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"fs"</code> - integer sample rate in samples/second 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"Nsamples"</code> - integer number of samples in the file (in each
channel) 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"wavetype"</code> - string with the name file type read. 
</p>

Other keys and values in the `info` dictionary are recorder-specific and
depend on the `wavetype` value.

### Recorder-Specific info keys and values

#### AudioMoth

Recordings identified as AudioMoth recordings have `info["wavetype"]` set
to "AudioMoth".  In addition to the standard `info` parameters, the following metadata parameters are added:

<p style="margin-left: 3em; text-indent: -2em;">
<code>"ICMT"</code> -  string with the contents of the ICMT subchunk.
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"IART"</code> -  string with the contents of the IART subchunk.
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"datestring"</code> -  string with the date and time of the
beginning of the recording in ISO 8601 format.
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"voltage"</code> -  string with the battery voltage at the
beginning of the recording.
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"gain"</code> -  string with the AudioMoth gain setting for
            recording.
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"serial number"</code> -  string with the serial number of the
AudioMoth recording device. 
</p>


#### Decimus

Recordings identified as Decimus recordings have `info["wavetype"]` set
to "decimus". Otherwise, `info` contains only the standard info keys and
values.

#### Generic

Recordings classified as generic have `info["wavetype"]` set to
"generic", and `info` contains only the standard info keys and values.

#### icListen

Recordings identified as icListen recordings have `info["wavetype"]` set
to "icListen". In addition, each key/value pair encoded in the INFO
chunk in the file is added to `info`. See the icListen documentation
for details on these parameters. 

The value `info["cal"]` contains a float64 calibration value for the data.
Multiply data samples by this value to obtain calibrated values in
micropascals.

#### Zoom

Recordings identified as Zoom recordings have `info["wavetype"]` set to
"zoom". The following information encoded in the bext chunk is added to
`info` as keys and values. (See Zoom documentation for details.)
    
<p style="margin-left: 3em; text-indent: -2em;">
<code>"CodingHistory"</code> - coding history string 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"desc"</code> - recording description string 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"LoudnessRange"</code> - int16 recording loudness range value 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"LoudnessValue"</code> - int16 recording loudness value 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"MaxMomentaryLoudness"</code> - int16 recording maximum momentary
loudness value 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"MaxShortTermLoudness"</code> - int16 recording maximum short term
loudness value
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"MaxTruePeakLevel"</code> - int16 recording maximum maximum true
peak level 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"OriginationDate"</code> - recording origination date string 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"OriginationTime"</code> - recording origination time string 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"Originator"</code> - recording originator string 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"OriginatorReference"</code> - recording originator reference
string 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"TimeReferenceHigh"</code> - int32 time of high sample in
recording 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"TimeReferenceLow"</code> - int32 time of low sample in recording 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"UMID"</code> - UMID string 
</p>

<p style="margin-left: 3em; text-indent: -2em;">
<code>"Version"</code> - int16 version number 
</p>

The contents in the entire iXML block are stored in `info["iXML"]` as a
string.

## levels Module

The module `levels` contains the functions `spl`, `sel`, `spl_wav`, `spl_wav_dir`, `spl_wav_files`, `A_weighting`, `M_weighting`, and `weight`. Each funcion contains a detailed usage message.

### spl

Return an array of sampled sound pressure levels using time constant T.  

#### Usage
    
        SPL = spl(data, fs, weighting='A', tconst=0.125, pref=20.0)
    
#### Input Parameters

`data`: an array of sampled sound pressures.
    
`fs`: sampling frequency in hertz.

`weighting`: type of weighting to use. This parameter can be a string to represent
preset values 'A' for A-weighting, 'M' for M-weighting (see
documentation on the function weight() to set frequency parameters). It
can also be a function that provides digital filter parameters to the
weight() function. For no weighting, use weighting = 1. The default is
'A' weighting.

    tconst: time constant. Defaullts to 0.125 s (fast). This parameter can be 
        the value in seconds or preset values given with the strings 'Fast' 
        (0.125 s), 'Slow' (1.000 s), or "Impulse' (0.035 s).
    pref: reference pressure. Defaults to 20.0 (micropascals, standard for
        atmospheric sounds). Use 1.0 (micropascals) for underwater sounds.
    cal: calibration factor of the recording. This is the value that
        converts data samples to the appropriate pressure units
        (micropascals). The default value is 1 (no calibration
        adjustment).
    pms: an initial value for the mean square pressure 'historical'
        value for time constant. Use this to continue the calculation
        from another recording. Defaults to 0.0.
    pms_return: whether or not to return the mean square pressure value for 
        subsequent calculations. Defaults to False.
    
    Output
    
    SPL: a numpy array of sampled sound pressure levels corresponding to the 
        the same sampling times a the elements of data.  Note that the initial 
        elements SPL[i] are based on a truncated history because they only use 
        pressure values from data[i] back to data[0].
    pms: The mean square sound pressure for use in subsequent
        calculations such as the recording continuing in another
        file. Only returned if the input parameter pms_return is True.

## Usage Example

Read data from a single-channel file and plot it vs. time.

    >>> from __future__ import division
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from SVSound import wavefile
    >>> info, data = wavefile.read('filename.wav')
    >>> info['chan']
    1
    >>> times = np.arange(data.size / info['fs'])
    >>> plt.plot(times, data)
    ...

Note that the data in a multichannel recording has rows for each channel, so `data[0]` is the first channel, `data[1]` the second channel, etc.
    
