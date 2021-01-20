# spec2nii
![PyPI](https://img.shields.io/pypi/v/spec2nii)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/spec2nii)  
A program for multi-format conversion of in vivo MRS to the [NIfTI MRS format](https://github.com/wexeee/mrs_nifti_standard).

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
| GE            | .7 (pfile)     | Yes | No  | No                    |
| UIH DICOM     | .dcm           | Yes | Yes | Yes                   |
| LCModel       | .RAW           | Yes | No  | No                    |
| jMRUI         | .txt           | Yes | No  | No                    |
| jMRUI         | .mrui          | Yes | No  | No                    |
| ASCII         | .txt           | Yes | No  | No                    |

## Instructions
spec2nii is called on the command line, and the conversion file type is specified with a subcommand.

The -f and -o options specify output file name and directory respectively for all formats.

If -j is specified the NIfTI MRS header extension will also be generated as a JSON side car file.

By default, spec2nii generates NIfTI files using the NIfTI-2 header format. Use the `--nifti1` option to generate files using the NIfTI-1 format.

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

### GE (limited support)
`spec2nii GE FILE`

### Philips (limited support)
`spec2nii GE SDAT_FILE SPAR_FILE`

### Text/LCModel/jMRUI
Conversion from processed formats.
All take an optional -a argument to specify a text file containing a 4x4 affine matrix specifying orientation information.

The text format requires additional information, namely imaging frequency in MHz and bandwidth in hertz.

`spec2nii raw -a AFFINE_FILE FILE`  
`spec2nii jmrui -a AFFINE_FILE FILE`  
`spec2nii text -a AFFINE_FILE -i imaging_freq -b bandwidth FILE`

