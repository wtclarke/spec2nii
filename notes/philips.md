# Philips data formats

## SPAR/SDAT

## .data/.list

## DICOM
There are two formats of Philips DICOM files. The two types can be identified as described by [Chris Rorden](https://github.com/rordenlab/dcm2niix/tree/master/Philips#non-image-dicoms)

The newer format `MR Spectroscopy Storage` can be identified by the `MediaStorageSOPClassUID (0002,0002)` being set to `1.2.840.10008.5.1.4.1.1.4.2`. The older format `Private MR Spectrum Storage` is set to `1.3.46.670589.11.0.0.12.1`.

There are many differences between the format of the two types, not least data storage location, header structure, and the use of multiple files versus one to store the complete data set in.

For DICOM tag interpretation in the old format see https://github.com/malaterre/dicom-private-dicts/blob/master/PMS-R32-dict.txt

### Header description for enhanced
 - (2005,1597) – tells if it’s an editing sequence
 - (2005,1598) – describes the editing type of the spectrum (0 indicates OFF, 1 indicates ON)
 - (2005,1304) – spectrum_mix_no (tells if its water data)
 - (2005,1357) – Sample Frequency (BW / spectral width)
 - (2001,1083) – Imaging Frequency (main field resonance frequency in MHz)
 - (2005,1310) – Spectrum Echo time
 - (0018,0080) – Repetition Time
 
 - (2005,1054) – volume_angulation_ap
 - (2005,1055) – volume_angulation_fh
 - (2005,1056) – volume_angulation_rl
 - (2005,105A) – volume_offcentre_ap
 - (2005,105B) – volume_offcentre_fh
 - (2005,105C) – volume_offcentre_rl
 
Under (2005,1085) you can find the FOV size
 - (2005,1057) – FOV in AP
 - (2005,1058) – FOV in FH
 - (2005,1059) – FOV in RL