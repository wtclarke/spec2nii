# Philips data formats

## SPAR/SDAT

## .data/.list

## DICOM
There are two formats of Philips DICOM files. The two types can be identified as described by [Chris Rorden](https://github.com/rordenlab/dcm2niix/tree/master/Philips#non-image-dicoms)

The newer format `MR Spectroscopy Storage` can be identified by the `MediaStorageSOPClassUID (0002,0002)` being set to `1.2.840.10008.5.1.4.1.1.4.2`. The older format `Private MR Spectrum Storage` is set to `1.3.46.670589.11.0.0.12.1`.

There are many differences between the format of the two types, not least data storage location, header structure, and the use of multiple files versus one to store the complete data set in.

For DICOM tag interpretation in the old format see https://github.com/malaterre/dicom-private-dicts/blob/master/PMS-R32-dict.txt 