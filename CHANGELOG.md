This document contains the Spec2nii release history in reverse chronological order.

0.5.0 (Friday 2nd December 2022)
--------------------------------
- NIfTI-MRS version 0.6
- Handle locale specific decimal separators in rda format.
- Added special-case handling of the mgs_svs_ed(_universal) for twix data inputs.

0.4.9 (Friday 4th November 2022)
--------------------------------
- Handle FID sequence in twix converter.
- Code lint updates
- Updated versioneer.

0.4.8 (Sunday 2nd October 2022)
-------------------------------
- Updated citation information with publication of NIfTI-MRS paper.
- Updated README information.
- Added check that file had been written successfully.
- Minor under the hood changes

0.4.7 (Sunday 11th September 2022)
----------------------------------
- Add tests for previous Siemens SOP UID DICOM types from VB line scanners.
- Fixed bug in Siemens RDA header read. Still no orientation tests for RDA format.
- Siemens VX line dicom files will now sum arrays of TE values.
- Added better ID of the Siemens implementation of the standardised SLASER sequence (as DICOM).
- Added provisional support for the GE slaser_cni sequence.
- Added handling of LCModel RAW files lacking header information.
- Included code and spelling fixes kindly contributed by Dimitri Papadopoulos

0.4.6 (Friday 8th July 2022)
----------------------------
- Fix GE p file incorrect TR and TE. Were expressed in microseconds.

0.4.5 (Friday 8th July 2022)
----------------------------
- Improved handling of Siemens DICOM files. Now explicitly handle Siemens Syngo Non Image Storage or MRSpectroscopyStorage DICOM types, regardless of scanner baseline.

0.4.4 (Thursday 28th April 2022)
--------------------------------
- Improved handling of Philips enhanced DICOM
- WIP handling of Philips classic DICOM, more test data required.
- SPAR/SDAT Philips data can now be reshaped ont eh command line using the -s and -t commands.
- Corrected phase convention for Siemens NumarisX (XA) line DICOM data. This was reversed (conjugated) compared to VX line.
