This document contains the Spec2nii release history in reverse chronological order.

0.7.3 (Tuesday 12th March 2024)
-------------------------------
- Siemens .rda format now had corrected and validated orientations (tested on VE11 baseline).
- Siemens .rda format now handles MRSI/CSI data and matches DICOM output. Validated on VE11 baseline data.
- Fixes in Siemens Twix special case for universal editing sequence (HERMES conditions).
- Added handling of custom Bruker sequences `mt_sLASER`, `mt_MEGA_sLASER_V35` and `cl_STELASER_PA360_b`.
- Philips vendor MEGA-PRESS handled through DICOM pathway. Thanks to Sandeep Ganji and Yansong Zhao for their help.

0.7.2 (Thursday 7th December 2023)
----------------------------------
- SpectralWidth now added to header extension automatically to match bids specification.
- NIfTI-MRS V0.8 now generated.
- Better handling of philips spar/sdat tags and singleton dimensions.
- Fixed bug where no name was set when a mixed folder of imaging and spectroscopy dicom was provided.

0.7.1 (Tuesday 7th November 2023)
---------------------------------
- The --anon flag can be passed with any call to anonymise after writing files.
- The Siemens enhanced dicom filetype pathway now handles CSI data.
- Fixed issue with RDA files having latin1 encoding. Thanks to gaunab on github. Fixes Issue #96.
- Now support GE p-files up to version 30.0.

0.7.0 (Saturday 5th August 2023)
--------------------------------
- Fixed a bug in Philips Classic DICOM orientations (supplementing the fixes to Enhanced DICOM in `0.6.11`)
- Updated Hyper sequence special case for Siemens twix and GE.

0.6.11 (Friday 28th July 2023)
------------------------------
- Fixed a bug in Philips spar/sdat orientations - this will affect voxels with rotations in more than one axes.
- Fixed a bug in Philips DICOM orientations - this will affect voxels with rotations in more than one axes.
- Python 3.7 now in [end of life](https://devguide.python.org/versions/) status and is no longer supported.
- `spec2nii dump/extract/insert` can now be used to inspect and fix non-compliant NIfTI-MRS files.

0.6.10 (Monday 10th July 2023)
------------------------------
- Stop tests being vendored as a top-level package. Bug introduced in `0.6.9`.

0.6.9 (Friday 7th July 2023)
----------------------------
- Add handling for the `fidall` sequence in GE.
- Add handling for the `slaser` sequence in GE.
- Update to nifti-mrs-tools 1.0.0 API.
- Minor code fixes contributed by the community.

0.6.8 (Wednesday 22nd March 2023)
---------------------------------
- Added handling for the GE jpress sequence (for JH HURCULES sequence)
- Added twix handling of the HYPER (smm_svs_herc_hyper) sequence (added by Aaron Gudmundson)
- Better handling of partial acquisitions of Siemens and Philips HYPER sequence.

0.6.7 (Wednesday 15th March 2023)
---------------------------------
- `spec2nii insert` can now insert a compliant header object into a non-compliant NIfTI file to generate a new NIfTI-MRS file.
- `spec2nii insert` can now update dwelltime using the optional `--dwelltime` argument.

0.6.6 (Thursday 9th March 2023)
-------------------------------
- Added in ability to generate empty NIfTI-MRS for twix-pathway MRSI scans.

0.6.5 (Wednesday 8th February 2023)
-----------------------------------
- Fixed bug in philips and rda metadata translation functions.

0.6.4 (Tuesday 7th February 2023)
---------------------------------
- Added first pass at new `spec2nii auto` feature with automatic conversion for some formats.

0.6.3 (Monday 6th February 2023)
--------------------------------
- Added handling for GE presscsi and probe-s sequences.

0.6.2 (Sunday 5th February 2023)
----------------------------------
- Handle HYPER references in SPAR/SDAT pipeline
- Handle HURCULES/HERMES (smm_svs_herc) sequence in XA twix format.
- Changed behaviour of Siemens DICOM `spec2nii dicom` to recursively glob directory argument.
- Better error output when encountering MR Image SOPClassUID.

0.6.1 (Wednesday 18th January 2023)
-----------------------------------
- Fixed conjugation issue introduced by new nifti-mrs package dependency
- SPAR/SDAT pipeline now handles HYPER special case.
- Data/list pipeline now handles HYPER special case.
- Fixed issue with XA Twix PatientSex and TxOffset attributes.
- Re-enable Bruker conversion.

0.6.0 (Wednesday 11th January 2023)
-----------------------------------
- NIfTI-MRS creation/handling/verification now performed by nifti-mrs package.

0.5.0 (Friday 2nd December 2022)
--------------------------------
- NIfTI-MRS version 0.6
- Handle locale specific decimal separators in rda format.
- Added special-case handling of the mgs_svs_ed(_universal) for twix data inputs.
- Modified repository structure. Siemens and Philips files now have dedicated directories.
- Support Python 3.11
- Enable automatic Pypi upload and move to pyproject.toml file build.
- Bruker conversion currently limited to numpy installations <1.20.0. Awaiting brukerapi package update.

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
