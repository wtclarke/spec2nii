"""Module containing other miscellaneous functions

Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2021 University of Oxford
"""
import json
import pprint

import nibabel as nib

from spec2nii import nifti_mrs


def dump_headers(args):
    """Print the nifti headers and header extensions

    :param args: Command line argumentss
    """
    # Load file
    nifti_mrs_img = nib.load(args.file)

    # Print NIfTI header
    print('NIfTI header: ')
    print(nifti_mrs_img.header)

    hdr_ext_codes = nifti_mrs_img.header.extensions.get_codes()
    hdr_ext = json.loads(nifti_mrs_img.header.extensions[hdr_ext_codes.index(44)].get_content())

    print('\nNIfTI-MRS header extension: ')
    pp = pprint.PrettyPrinter()
    pp.pprint(hdr_ext)


def extract_hdr_ext(args):
    """Generate a json sidecar file for input file

    :param args: Command line argumentss
    """
    nifti_mrs_img = nib.load(args.file)

    hdr_ext_codes = nifti_mrs_img.header.extensions.get_codes()

    if args.outdir:
        args.outdir.mkdir(exist_ok=True, parents=True)
        out_json = args.outdir / (args.file.with_suffix('').with_suffix('').name + '.json')
    else:
        out_json = args.file.parent / (args.file.with_suffix('').with_suffix('').name + '.json')

    with open(out_json, 'w') as fp:
        json.dump(json.loads(nifti_mrs_img.header.extensions[hdr_ext_codes.index(44)].get_content()), fp, indent=4)


def insert_hdr_ext(args):
    """Function for anonymising input NIfTI-MRS files.

    :param args: Command line arguments parsed in spec2nii.py
    :return: List of anonymised images
    :rtype: [nib.nifti2.Nifti2Image,]
    :return: List of output names
    :rtype: [str,]
    """
    # Load data
    nifti_mrs_img = nib.load(args.file)

    with open(args.json_file) as jf:
        new_hdr = json.load(jf)

    # Make new NIfTI-MRS image
    mod_out = nifti_mrs.NIfTI_MRS(nifti_mrs_img.get_fdata(dtype=nifti_mrs_img.get_data_dtype()),
                                  nifti_mrs_img.affine,
                                  nifti_mrs_img.header['pixdim'][4],
                                  new_hdr)

    # Process output name.
    if args.fileout:
        fname_out = [args.fileout, ]
    else:
        fname_out = [args.file.with_suffix('').with_suffix('').name, ]

    return [mod_out, ], fname_out
