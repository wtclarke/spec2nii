# spec2nii
![PyPI](https://img.shields.io/pypi/v/spec2nii)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/spec2nii)  
A program for multi-format in vivo MR spectroscopy conversion to NIFTI.

## Installation
`pip install spec2nii`

## Currently supported formats
This table lists the currently supported formats. I have very limited experience with Philips and GE formats. Please get in touch if you are willing to help add to this list and/or supply validation data.

| Format        | File extension | SVS | CSI | Automatic orientation |
|---------------|----------------|-----|-----|-----------------------|
| Siemens Twix  | .dat           | Yes | No  | Yes                   |
| Siemens DICOM | .ima / .dcm    | Yes | Yes | Yes                   |
| Philips       | .SPAR/.SDAT    | Yes | No  | No                    |
| GE            | .7 (pfile)     | Yes | No  | No                    |
| UIH DICOM     | .dcm           | Yes | No  | Yes                   |
| LCModel       | .RAW           | Yes | No  | No                    |
| jMRUI         | .txt           | Yes | No  | No                    |
| ASCII         | .txt           | Yes | No  | No                    |

## Instructions
spec2nii is called on the comandline, and the conversion file type is specified with a subcommand.

The -f and -o options specify output file name and directory respectively for all formats.

If -j is specified a sidecar JSON file containing usefull meta data will be generated from the file headers. This is similar to the BIDs format.

### Twix
`spec2nii twix`
Call `spec2nii twix -v FILE` to view a list of contained MDH flags. -m can be used to specify which multi-raid file to convert if used on VE data.

Call with the -e flag to specify which MDH flag to convert. e.g.  
`spec2nii twix -e image FILE`

### DICOM
`spec2nii dicom DCM_FILE_or_DIR`

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

