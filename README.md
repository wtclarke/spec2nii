# spec2nii
![PyPI](https://img.shields.io/pypi/v/spec2nii)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/spec2nii)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5907960.svg)](https://doi.org/10.5281/zenodo.5907960)

A program for multi-format conversion of in vivo MRS to the [NIfTI-MRS format](https://github.com/wexeee/mrs_nifti_standard).

## About

This program was inspired by the imaging DICOM to NIfTI converter [dcm2niix](https://github.com/rordenlab/dcm2niix) developed by Chris Rorden. All MRS(I) orientations are tested with images converted using dcm2niix. I consider the combination of images converted using dcm2niix and displayed in [FSLeyes](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLeyes) the de facto standard.

## Citing Spec2nii

If you use `spec2nii` or the [NIfTI-MRS](https://github.com/wexeee/mrs_nifti_standard) format in your work please cite:

    `Clarke WT, Bell TK, Emir UE, Mikkelsen M, Oeltzschner G, Shamaei A, Soher BJ, Wilson M. NIfTI-MRS: A standard data format for magnetic resonance spectroscopy. Magn Reson Med. 2022. doi: 10.1002/mrm.29418.`

## Visualising output

Visualisation of MRS data converted by spec2nii to NIfTI-MRS can be carried out with a recent (>0.31.0) version of FSLeyes. A FSLeyes plugin for NIfTI-MRS is now available:
- [Gitlab](https://git.fmrib.ox.ac.uk/wclarke/fsleyes-plugin-mrs),
- `conda install -c conda-forge fsleyes-plugin-mrs`,
- `pip install fsleyes-plugin-mrs`.

## Installation
`conda install -c conda-forge spec2nii`
or
`pip install spec2nii`

### Installing Conda (option #1)
Miniconda can be installed by following the instructions on the [Conda website](https://docs.conda.io/en/latest/miniconda.html). To create a suitable environment run the following three commands after installing Conda.

```
    conda create -c conda-forge -n my_env python=3.8
    conda activate my_env
    conda install -c conda-forge spec2nii
```

## Currently supported formats
This table lists the currently supported formats. I have very limited experience with Philips and GE formats. Please get in touch if you are willing to help add to this list and/or supply validation data.

| Format        | File extension | SVS | MRSI| Automatic orientation |
|---------------|----------------|-----|-----|-----------------------|
| Siemens Twix  | .dat           | Yes | †   | Yes                   |
| Siemens DICOM | .ima / .dcm    | Yes | Yes | Yes                   |
| Siemens RDA   | .rda           | Yes | Yes | Yes                   |
| Philips       | .SPAR/.SDAT    | Yes | No  | Yes                   |
| Philips       | .data/.list    | Yes | No  | Yes                   |
| Philips DICOM | .dcm           | Yes | No  | Yes                   |
| GE            | .7 (pfile)     | Yes | Yes | Yes  (WIP)            |
| UIH DICOM     | .dcm           | Yes | Yes | Yes                   |
| Bruker        | 2dseq          | Yes | Yes | Yes                   |
| Bruker        | fid            | Yes | Yes | Yes (WIP)             |
| Varian        | fid            | Yes | No  | No (WIP)              |
| LCModel       | .RAW           | Yes | No  | No                    |
| jMRUI         | .txt           | Yes | No  | No                    |
| jMRUI         | .mrui          | Yes | No  | No                    |
| ASCII         | .txt           | Yes | No  | No                    |

† Partial handling - see section on Twix pathway for MRSI handling.

## Instructions
spec2nii is called on the command line, and the conversion file type is specified with a subcommand.

The -f and -o options specify output file name and directory respectively for all formats.

If -j is specified the NIfTI MRS header extension will also be generated as a JSON side car file.

By default, spec2nii generates NIfTI files using the NIfTI-2 header format. Use the `--nifti1` option to generate files using the NIfTI-1 format.

Manual overrides can be provided for incorrectly interpreted required header fields, namely SpectrometerFrequency, ResonantNucleus and dwell-time, by using the `--override_frequency`, `--override_nucleus`, and `--override_dwelltime` command line options.

### Automatic detection
`spec2nii auto FILE` will attempt an automatic conversion of the following formats: Twix, RDA, SPAR/SDAT, GE p-file, DICOM. Note that many features of the individual converters are not implemented in this automatic pathway. This feature should be regarded as somewhat experimental. For finer-grained control see the specific subcommands listed below.

### Twix
Call `spec2nii twix -v FILE` to view a list of contained MDH flags. -m can be used to specify which multi-raid file to convert if used on VE data.

Call with the -e flag to specify which MDH flag to convert. e.g.
`spec2nii twix -e image FILE`

Twix format loop variables (e.g. `Ave` or `ida`) can be assigned to specific NIfTI dimensions using the `-d{5,6,7}` command line options. NIfTI MRS dimension tags (e.g. `DIM_COIL`) can be specified using the `-t{5,6,7}` command line options.

As `spec2nii` is __not__ a reconstruction program, it cannot convert MRSI data. Far too little information is held in the twix headers to reconstruct arbitrary k,t-space data. However, if passed a file containing MRSI data `spec2nii` will attempt to create an empty NIfTI-MRS file with the correct orientation information, data shape, and header information. This empty file can then have data inserted from an offline reconstruction routine.

### Siemens DICOM
`spec2nii dicom DCM_FILE_or_DIR`

NIfTI MRS dimension tags (e.g. `DIM_COIL`) can be specified using the `-t` command line argument.

### Siemens RDA
`spec2nii rda RDA_FILE`

Compatible with CSI and SVS data. Validated to be the same data and orientation information as DICOM output on VE baselines.

### UIH DICOM
`spec2nii uih DCM_FILE_or_DIR`

NIfTI MRS dimension tags (e.g. `DIM_COIL`) can be specified using the `-t` command line argument.

### GE
`spec2nii ge FILE`

### Philips (SPAR/SDAT)
`spec2nii philips SDAT_FILE SPAR_FILE`

Two optional arguments are available for the SPAR/SDAT pathway:
- `-t/--tags` allows the user to specify the dimension tags for each of the higher dimensions (up to three).
- `-s/--shape` allows the user to perform Numpy style reshaping of multiple transients. By default (without specifying a shape) all transients will be listed in a single 5th dimension.

### Philips (data/list)
Must be provided along side a matching SPAR file.
`spec2nii philips_dl DATA_FILE LIST_FILE SPAR_FILE`

### Philips DICOM
`spec2nii philips_dcm DCM_FILE_or_DIR`

NIfTI MRS dimension tags (e.g. `DIM_COIL`) can be specified using the `-t` command line argument.

Generates separate reference file if present.

Both classic and enhanced DICOM is handled. Well tested on the vendors' own PRESS and MEGA-PRESS sequence.

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

Note that the varian file format is very flexible -- the binary `fid` itself essentially is a long 2D list of (`complex_points * everything_else`), and the current code makes several significant assumptions about how that should be interpreted and reshaped.
In particular, if you are using a sequence derived from something different to either `spuls`, `s2pul`, `press`, or `steam`, it is quite likely that this will not work. Edit `varian_importer.py` and add cases based on your `seqfil` as appropriate.
It is assumed that the `comment` parameter should be the patient's name.

(Further bells and whistles pending; Written by Jack J. Miller <jack.miller@physics.org>)


### Text/LCModel/jMRUI
Conversion from processed formats.
All take an optional -a argument to specify a text file containing a 4x4 affine matrix specifying orientation information.

The text format requires additional information, namely imaging frequency in MHz and bandwidth in hertz.

`spec2nii raw -a AFFINE_FILE FILE`
`spec2nii jmrui -a AFFINE_FILE FILE`
`spec2nii text -a AFFINE_FILE -i imaging_freq -b bandwidth FILE`

### Other functions
Anonymise the NIfTI-MRS file. All standard-defined keys marked as sensitive will be removed. User defined parameters marked as `private_` will also be removed. Use the `-v` flag to view the removed header keys. The `-r` argument may be used (repeatedly) to remove additional keys from the header extension manually.
`spec2nii anon FILE`

Dump the NIfTI headers and header extension to Stdout.
`spec2nii dump FILE`

Produce a json file containing the header extension as a separate file from a NIfTI-MRS file.
`spec2nii extract FILE`

Overwrite the header extension in a NIfTI-MRS file using a separate json formatted file.
`spec2nii insert FILE JSON_FILE`

## Contributors & contributing
This program was written by Will Clarke, University of Oxford. Contributions to add new file formats or improve existing ones are very welcome. Please fork the repository and request changes using a merge (pull) request. I ask that test data and tests are included with any submission.

Particular thanks go to Tomáš Pšorn for contributing the Bruker interface, and to Jack Miller for the Varian interface.

Elements of the Varian reader come from [NMR glue](https://github.com/jjhelmus/nmrglue/), if you use the varian components in your research please cite J.J. Helmus, C.P. Jaroniec, Nmrglue: An open source Python package for the
analysis of multidimensional NMR data, J. Biomol. NMR 2013, 55, 355-367. doi: 10.1007/s10858-013-9718-x

Some GE test data comes from the [BIG GABA](https://www.nitrc.org/projects/biggaba/) dataset which was funded by NIH grant R01 EB016089. Please see Mikkelsen M et al. Big GABA: Edited MR spectroscopy at 24 research sites. NeuroImage 2017;159:32–45. doi: 10.1016/j.neuroimage.2017.07.021 for more information.
