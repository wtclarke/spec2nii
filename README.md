# spec2nii
![PyPI](https://img.shields.io/pypi/v/spec2nii)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/spec2nii)  
A program for multi-format conversion of in vivo MRS to the [NIfTI-MRS format](https://github.com/wexeee/mrs_nifti_standard).  

This program was inspired by the imaging DICOM to NIfTI converter [dcm2niix](https://github.com/rordenlab/dcm2niix) developed by Chris Rorden. All MRS(I) orientations are tested with images converted using dcm2niix. I consider the combination of images converted using dcm2niix and displayed in [FSLeyes](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLeyes) the de facto standard.

Visualisation of MRS converted by spec2nii can be carried out with a recent (>0.31.0) version of FSLeyes. A FSLeyes plugin for NIfTI-MRS is in development.
## Installation
`conda install -c conda-forge spec2nii`  
or  
`pip install spec2nii`

## Currently supported formats
This table lists the currently supported formats. I have very limited experience with Philips and GE formats. Please get in touch if you are willing to help add to this list and/or supply validation data.

| Format        | File extension | SVS | CSI | Automatic orientation |
|---------------|----------------|-----|-----|-----------------------|
| Siemens Twix  | .dat           | Yes | No  | Yes                   |
| Siemens DICOM | .ima / .dcm    | Yes | Yes | Yes                   |
| Philips       | .SPAR/.SDAT    | Yes | No  | Yes                   |
| Philips       | .data/.list    | Yes | No  | Yes                   |
| Philips DICOM | .dcm           | Yes | No  | Yes (WIP)             |
| GE            | .7 (pfile)     | Yes | Yes | Yes                   |
| UIH DICOM     | .dcm           | Yes | Yes | Yes                   |
| Bruker        | 2dseq          | Yes | Yes | Yes                   |
| Bruker        | fid            | Yes | Yes | Yes (WIP)             |
| Varian        | fid            | Yes | No  | No (WIP)              |
| LCModel       | .RAW           | Yes | No  | No                    |
| jMRUI         | .txt           | Yes | No  | No                    |
| jMRUI         | .mrui          | Yes | No  | No                    |
| ASCII         | .txt           | Yes | No  | No                    |

## Instructions
spec2nii is called on the command line, and the conversion file type is specified with a subcommand.

The -f and -o options specify output file name and directory respectively for all formats.

If -j is specified the NIfTI MRS header extension will also be generated as a JSON side car file.

By default, spec2nii generates NIfTI files using the NIfTI-2 header format. Use the `--nifti1` option to generate files using the NIfTI-1 format.

Manual overrides can be provided for incorrectly interpreted required header fields, namely SpectrometerFrequency, ResonantNucleus and dwell-time, by using the `--override_frequency`, `--override_nucleus`, and `--override_dwelltime` command line options.

### Twix
Call `spec2nii twix -v FILE` to view a list of contained MDH flags. -m can be used to specify which multi-raid file to convert if used on VE data.

Call with the -e flag to specify which MDH flag to convert. e.g.  
`spec2nii twix -e image FILE`

Twix format loop variables (e.g. `Ave` or `ida`) can be assigned to specific NIfTI dimensions using the `-d{5,6,7}` command line options. NIfTI MRS dimension tags (e.g. `DIM_COIL`) can be specified using the `-t{5,6,7}` command line options.

### Siemens DICOM
`spec2nii dicom DCM_FILE_or_DIR`

NIfTI MRS dimension tags (e.g. `DIM_COIL`) can be specified using the `-t` command line argument.

### UIH DICOM
`spec2nii uih DCM_FILE_or_DIR`

NIfTI MRS dimension tags (e.g. `DIM_COIL`) can be specified using the `-t` command line argument.

### GE
`spec2nii ge FILE`

### Philips (SPAR/SDAT)
`spec2nii philips SDAT_FILE SPAR_FILE`

### Philips (data/list)
Must be provided along side a matching SPAR file.  
`spec2nii philips_dl DATA_FILE LIST_FILE SPAR_FILE`

### Philips DICOM
`spec2nii philips_dcm DCM_FILE_or_DIR`

NIfTI MRS dimension tags (e.g. `DIM_COIL`) can be specified using the `-t` command line argument.

Generates separate reference file.

### Bruker (2dseq/fid)
`spec2nii bruker -m 2DSEQ 2DSEQ_FILE_or_DIR`  
`spec2nii bruker -m FID FID_FILE_or_DIR`

Use the `-d` option to dump the header files (method and acqp for fid, visu_pars for 2dseq) into the header extension.

Additional filters can be added by defining additional queries using the `-q` flag.

Bruker conversion is powered by the [BrukerAPI package](https://github.com/isi-nmr/brukerapi-python) written by Tomas Psorn.

### Varian 
`spec2nii varian /path/to/fid.fid`  
where fid.fid is a Varian fid directory containing a fid and procpar file.  
Use the `-d` option to dump the procpar header file contents into the header extension.  
Use the `-t` option to set an alternative dimension tag for the 6th dimension (default = `DIM_DYN`).  

(Bells and whistles pending -- this only really works with 1D spectra that may change over time and may be received on multiple coils)  
Written by Jack J. Miller (jack.miller@physics.org) 


### Text/LCModel/jMRUI
Conversion from processed formats.
All take an optional -a argument to specify a text file containing a 4x4 affine matrix specifying orientation information.

The text format requires additional information, namely imaging frequency in MHz and bandwidth in hertz.

`spec2nii raw -a AFFINE_FILE FILE`  
`spec2nii jmrui -a AFFINE_FILE FILE`  
`spec2nii text -a AFFINE_FILE -i imaging_freq -b bandwidth FILE`

