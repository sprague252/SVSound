"""This subpackage contains modules with format specifications for several proprietary Broadcast Wave file formats. Currently supported formats include:

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
__all__ = ["decimus", "generic", "icListen", "zoom"]

