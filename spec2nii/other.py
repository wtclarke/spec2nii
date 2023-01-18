"""Module containing other miscellaneous functions

Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2021 University of Oxford
"""
import json
import pprint

from nifti_mrs.nifti_mrs import NIFTI_MRS


def dump_headers(args):
    """Print the nifti headers and header extensions

    :param args: Command line arguments
    """
    # Load file
    nifti_mrs_img = NIFTI_MRS(args.file)

    # Print NIfTI header
    print('NIfTI header: ')
    print(nifti_mrs_img.header)

    print('\nNIfTI-MRS header extension: ')
    pp = pprint.PrettyPrinter()
    pp.pprint(nifti_mrs_img.hdr_ext.to_dict())


def extract_hdr_ext(args):
    """Generate a json sidecar file for input file

    :param args: Command line arguments
    """
    nifti_mrs_img = NIFTI_MRS(args.file)

    if args.outdir:
        args.outdir.mkdir(exist_ok=True, parents=True)
        out_json = args.outdir / (args.file.with_suffix('').with_suffix('').name + '.json')
    else:
        out_json = args.file.parent / (args.file.with_suffix('').with_suffix('').name + '.json')

    with open(out_json, 'w') as fp:
        json.dump(nifti_mrs_img.hdr_ext.to_dict(), fp, indent=4)


def insert_hdr_ext(args):
    """Function for inserting a new header into NIfTI-MRS files.

    :param args: Command line arguments parsed in spec2nii.py
    :return: Modified NIfTI-MRS file
    :rtype: [nifti_mrs.nifti_mrs.NIFTI_MRS,]
    :return: List of output names
    :rtype: [str,]
    """
    # Load data
    nifti_mrs_img = NIFTI_MRS(args.file)

    with open(args.json_file) as jf:
        new_hdr = json.load(jf)

    nifti_mrs_img.hdr_ext = new_hdr

    # Process output name.
    if args.fileout:
        fname_out = [args.fileout, ]
    else:
        fname_out = [args.file.with_suffix('').with_suffix('').name, ]

    return [nifti_mrs_img, ], fname_out
