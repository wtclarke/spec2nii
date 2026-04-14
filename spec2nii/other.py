"""Module containing other miscellaneous functions

Author:  William Clarke         <william.clarke@ndcn.ox.ac.uk>
         Vasilis Karlaftis      <vasilis.karlaftis@ndcn.ox.ac.uk>

Copyright (C) 2026 University of Oxford
"""
import json
import pprint

from nifti_mrs.nifti_mrs import NIFTI_MRS, NotNIFTI_MRS


def dump_headers(args):
    """Print the nifti headers and header extensions

    :param args: Command line arguments
    """
    # Load file
    nifti_mrs_img = NIFTI_MRS(args.file, validate_on_creation=False)

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
    nifti_mrs_img = NIFTI_MRS(args.file, validate_on_creation=False)

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

    with open(args.json_file) as jf:
        new_hdr = json.load(jf)

    # Load data
    try:
        nifti_mrs_img = NIFTI_MRS(args.file, validate_on_creation=False)
        nifti_mrs_img.hdr_ext = new_hdr

    except NotNIFTI_MRS:
        print(f'{args.file} is not compliant with the NIfTI-MRS standard, attempting to convert')
        from fsl.data.image import Image
        from nifti_mrs.create_nmrs import gen_nifti_mrs_hdr_ext
        from nifti_mrs.hdr_ext import Hdr_Ext

        nimg = Image(args.file)
        hdr_ext = Hdr_Ext.from_header_ext(new_hdr)
        nifti_mrs_img = gen_nifti_mrs_hdr_ext(
            nimg[:],
            nimg.header['pixdim'][4],
            hdr_ext,
            affine=nimg.getAffine('voxel', 'world'))

    if args.dwelltime is not None:
        nifti_mrs_img.dwelltime = args.dwelltime

    # Process output name.
    if args.fileout:
        fname_out = [args.fileout, ]
    else:
        fname_out = [args.file.with_suffix('').with_suffix('').name, ]

    return [nifti_mrs_img, ], fname_out


def clean_hdr_ext(args):
    """Function for updating invalid fields in NIfTI-MRS header.

    :param args: Command line arguments parsed in spec2nii.py
    :return: Modified NIfTI-MRS file
    :rtype: [nifti_mrs.nifti_mrs.NIFTI_MRS,]
    :return: List of output names
    :rtype: [str,]
    """
    import re
    from numpy import isclose, ndarray
    from nifti_mrs.definitions import standard_defined

    # Try reading the file as NIfTI-MRS
    try:
        nifti_mrs_img = NIFTI_MRS(args.file, validate_on_creation=False)
        print("'intent_name' is valid.")
    except NotNIFTI_MRS as exc:
        if str(exc) != 'NIFTI-MRS intent code not set.':
            raise
        # if failed, read as a regular NIfTI, update intent and re-save
        from fsl.data.image import Image
        from nifti_mrs.hdr_ext import Hdr_Ext
        img = Image(args.file)
        hdr_ext_codes = img.header.extensions.get_codes()
        if 44 not in hdr_ext_codes:
            raise NotNIFTI_MRS('NIFTI-MRS must have a header extension.')
        img._hdr_ext = Hdr_Ext.from_header_ext(
            img.header.extensions[hdr_ext_codes.index(44)].json())

        # 1. check that intent is of a valid format
        intent_ptrn = re.compile(r'mrs_v\d+_\d+')
        intent_str = img.header.get_intent()[2]
        if intent_str is None or intent_str == '' or intent_ptrn.match(intent_str) is None:
            import json
            from importlib.resources import files
            data_text = files('nifti_mrs.standard').joinpath('definitions.json').read_text(encoding='utf-8')
            json_def = json.loads(data_text)

            v_major = json_def['nifti_mrs_version']['major']
            v_minor = json_def['nifti_mrs_version']['minor']
            img.header['intent_name'] = f'mrs_v{v_major}_{v_minor}'.encode()
            print(f"'intent_name' is updated to 'mrs_v{v_major}_{v_minor}'.")

        # save file to output and re-read it
        args.outdir.mkdir(parents=True, exist_ok=True)
        if args.fileout:
            args.file = args.outdir / args.fileout
        else:
            args.file = args.outdir / args.file.with_suffix('').with_suffix('').name
        img.save(args.file)
        nifti_mrs_img = NIFTI_MRS(args.file, validate_on_creation=False)

    new_hdr = nifti_mrs_img.hdr_ext.copy()

    # 2. check that 'SpectrometerFrequency' is an array of floats
    val = new_hdr.SpectrometerFrequency
    if not isinstance(val, (list, tuple, ndarray)) or (isinstance(val, ndarray) and val.ndim == 0):
        val = [val]
    new_hdr.SpectrometerFrequency = [float(v) for v in val]
    print("'SpectrometerFrequency' format is updated.")

    # 3. check that 'ResonantNucleus' is an array of strings
    val = new_hdr.ResonantNucleus
    if not isinstance(val, (list, tuple, ndarray)) or (isinstance(val, ndarray) and val.ndim == 0):
        val = [val]
    new_hdr.ResonantNucleus = [str(v) for v in val]
    print("'ResonantNucleus' format is updated.")

    # 4. check that 'SpectralWidth' is equal to 1/pixdim[4]
    if 'SpectralWidth' in new_hdr:
        spec_width = new_hdr['SpectralWidth']
        if not isclose(spec_width, 1 / nifti_mrs_img.dwelltime, atol=1E-2):
            print("Warning: 'SpectralWidth' does not match '1/dwelltime'! "
                  f"Replacing with the latter: {1 / nifti_mrs_img.dwelltime}")
    else:
        print("'SpectralWidth' is added.")
    new_hdr.set_standard_def('SpectralWidth', 1 / nifti_mrs_img.dwelltime)

    # 5. check that user-defined fields are dictionary with a 'Description' field
    dim_re = re.compile(r"^dim_[567](_((info)|(header)))?$")
    for key in new_hdr:
        if key not in standard_defined\
           and key != "SpectrometerFrequency"\
           and key != "ResonantNucleus"\
           and not dim_re.match(key):
            # Must be user-defined, convert it to a dictionary with empty 'Description'
            if not isinstance(new_hdr[key], dict):
                val = new_hdr[key]
                new_hdr.set_user_def(key, val, '')
            elif "Description" not in new_hdr[key]:
                new_hdr[key]["Description"] = ''
            print(f"'{key}' is updated with empty 'Description'.")

    # update NIfTI-MRS header extension
    nifti_mrs_img.hdr_ext = new_hdr

    if args.fileout:
        fname_out = [args.fileout, ]
    else:
        fname_out = [args.file.with_suffix('').with_suffix('').name, ]

    return [nifti_mrs_img, ], fname_out
