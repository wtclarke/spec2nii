"""Function to anonymise existing NIfTI-MRS files.

Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2021 University of Oxford
"""

import re
import pprint

from nifti_mrs.definitions import standard_defined
from nifti_mrs.create_nmrs import gen_nifti_mrs_hdr_ext
from nifti_mrs.nifti_mrs import NIFTI_MRS
from nifti_mrs.hdr_ext import Hdr_Ext


def anon_nifti_mrs(nifti_mrs_img, extra_keys=None, verbose=False):
    """Function for anonymising an existing nifti-mrs object

    :param nifti_mrs_img: NIfTI-MRS object to anonymise
    :type nifti_mrs_img: nifti_mrs.nifti_mrs.NIFTI_MRS
    :param extra_keys: List of keys to explicitly remove, defaults to None
    :type extra_keys: List, optional
    :param verbose: Print verbose output, defaults to False
    :type verbose: bool, optional
    :return: Anonymised NIfTI-MRS object
    :rtype: nifti_mrs.nifti_mrs.NIFTI_MRS
    """

    # Extract header extension
    hdr_ext = nifti_mrs_img.hdr_ext.to_dict()

    # Loop through fields. Remove those which are marked for anonymisation
    # Either those fields which are standard-defined and are marked for anonymisation,
    # or user-defined which are marked 'private'.
    pvt_re = re.compile(r'private_')

    def iter_json(in_dict):
        out_dict = {}
        removed = {}
        for key, value in in_dict.items():
            # Explicitly set on command line
            if extra_keys\
                    and key in extra_keys:
                removed.update({key: value})
            # Standard defined
            elif key in standard_defined:
                if standard_defined[key][3]:
                    removed.update({key: value})
                else:
                    out_dict.update({key: value})
            # User defined
            else:
                # If marked as private at this level
                if pvt_re.match(key):
                    removed.update({key: value})
                # If nested dict
                elif isinstance(value, dict):
                    tmp_out, tmp_rem = iter_json(value)
                    if len(tmp_rem) > 0:
                        removed.update({key: tmp_rem})
                    out_dict.update({key: tmp_out})
                # Catch all other cases
                else:
                    out_dict.update({key: value})
        return out_dict, removed

    anon_hdr, removed_dict = iter_json(hdr_ext)

    if verbose:
        pp = pprint.PrettyPrinter(indent=4)
        print('\nThe following keys were removed:')
        pp.pprint(removed_dict)
        print('\nThe following keys were kept:')
        pp.pprint(anon_hdr)

    # Make new NIfTI-MRS image
    return gen_nifti_mrs_hdr_ext(
        nifti_mrs_img[:],
        nifti_mrs_img.dwelltime,
        Hdr_Ext.from_header_ext(anon_hdr),
        affine=nifti_mrs_img.getAffine('voxel', 'world'))


def anon_file(args):
    """Function for anonymising NIfTI-MRS files stored on drive.

    :param args: Command line arguments parsed in spec2nii.py
    :return: List of anonymised images
    :rtype: [nib.nifti2.Nifti2Image,]
    :return: List of output names
    :rtype: [str,]
    """

    # Process output name.
    if args.fileout:
        fname_out = [args.fileout, ]
    else:
        fname_out = [args.file.with_suffix('').with_suffix('').name, ]

    return [
        anon_nifti_mrs(
            NIFTI_MRS(args.file),
            verbose=args.verbose,
            extra_keys=args.remove), ], fname_out
