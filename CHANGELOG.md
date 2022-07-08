This document contains the FSL-MRS release history in reverse chronological order.

0.4.6 (Friday 8th July 2022)
----------------------------
- Fix GE p file incorrect TR and TE. Were expressed in microseconds.


0.4.5 (Friday 8th July 2022)
----------------------------
- Improved handling of Siemens DICOM files. Now explicitly handle Siemens Syngo Non Image Storage or MRSpectroscopyStorage DICOM types, irregardless of scanner baseline.

0.4.4 (Thursday 28th April 2022)
--------------------------------
- Improved handling of Philips enhanced DICOM
- WIP handling of Philips classic DICOM, more test data required.
- SPAR/SDAT Philips data can now be reshaped ont eh command line using the -s and -t commands.
- Corrected phase convention for Siemens NumarisX (XA) line DICOM data. This was reversed (conjugated) compared to VX line.
