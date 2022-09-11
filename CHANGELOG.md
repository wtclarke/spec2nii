This document contains the FSL-MRS release history in reverse chronological order.

0.4.7 (Sunday 11th September)
-----------------------------
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
