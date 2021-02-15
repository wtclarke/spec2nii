""" spec2nii - tool for conversion of various MRS data formats to NIFTI format.

spec2nii converts the following formats to NIFTI files.
Supporting SVS:
    Siemens "Twix" .dat format
    Siemens DICOM
    Philips SPAR/SDAT files
    GE p-files
    UIH DICOM
    LCModel RAW
    jMRUI text
    Plain text

Supporting CSI/MRSI:
    Siemens DICOM
    UIH DICOM

This module contains the main class to be called as a script (through the main function).

Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford
"""

import argparse
import sys
import os.path as op
from pathlib import Path
import json
from spec2nii import nifti_mrs
# There are case specific imports below


class spec2nii:
    def __init__(self):
        parser = argparse.ArgumentParser(description='Convert raw spectroscopy data to NIfTI format.')

        subparsers = parser.add_subparsers(title='subcommands',
                                           description='File types supported')

        def add_common_parameters(subparser):
            # Add options that are common to all subcommands
            subparser.add_argument('-j', '--json', help='Create json sidecar.', action='store_true')
            subparser.add_argument("-f", "--fileout", type=str,
                                   help="Output file base name (default = input file name)")
            subparser.add_argument("-o", "--outdir", type=Path,
                                   help="Output location (default = .)", default='.')
            subparser.add_argument('--nifti1', action='store_true')
            subparser.add_argument('--verbose', action='store_true')
            return subparser

        # Handle twix subcommand
        parser_twix = subparsers.add_parser('twix', help='Convert from Siemens .dat twix format.')
        parser_twix.add_argument('file', help='file to convert', type=Path)
        group = parser_twix.add_mutually_exclusive_group(required=True)
        group.add_argument("-v", "--view", help="View contents of twix file, no files converted", action='store_true')
        group.add_argument('-e', '--evalinfo', type=str, help='evalInfo flag to convert')
        parser_twix.add_argument("-m", "--multiraid", type=int, help="Select multiraid file to load"
                                 " (default = 2 i.e. 2nd file)", default=2)
        parser_twix.add_argument("-q", "--quiet", help="Suppress text output", action='store_true')
        for idx in range(5, 8):
            parser_twix.add_argument(f"-d{idx}", f"--dim{idx}", type=str, help=f"Specify dim {idx} loop counter.")
            parser_twix.add_argument(f"-t{idx}", f"--tag{idx}", type=str, help=f"Specify dim {idx} NIfTI MRS tag.")
        parser_twix = add_common_parameters(parser_twix)
        parser_twix.set_defaults(func=self.twix)

        # Handle dicom subcommand
        parser_dicom = subparsers.add_parser('dicom', help='Convert from Siemens DICOM format.')
        parser_dicom.add_argument('file', help='file or directory to convert', type=str)
        parser_dicom.add_argument("-t", "--tag", type=str, help="Specify NIfTI MRS tag used for 5th "
                                                                "dimension if multiple files are passed.")
        parser_dicom = add_common_parameters(parser_dicom)
        parser_dicom.set_defaults(func=self.dicom)

        # Handle UIH DICOM subcommand
        parser_uih_dicom = subparsers.add_parser('uih', help='Convert from UIH DICOM format.')
        parser_uih_dicom.add_argument('file', help='file or directory to convert', type=str)
        parser_uih_dicom.add_argument("-t", "--tag", type=str, help="Specify NIfTI MRS tag used for 5th "
                                                                    "dimension if multiple files are passed.")
        parser_uih_dicom = add_common_parameters(parser_uih_dicom)
        parser_uih_dicom.set_defaults(func=self.uih_dicom)

        # Handle philips subcommand
        parser_philips = subparsers.add_parser('philips', help='Convert from Philips spar/sdat format.')
        parser_philips.add_argument('sdat', help='SDAT file', type=Path)
        parser_philips.add_argument('spar', help='SPAR file', type=Path)
        parser_philips.add_argument("-t", "--tag", type=str, help="Specify NIfTI MRS tag used for 5th "
                                                                  "dimension if multiple transients are present.")
        parser_philips = add_common_parameters(parser_philips)
        parser_philips.set_defaults(func=self.philips)

        # Handle philips data/list subcommand
        parser_p_dl = subparsers.add_parser('philips_dl', help='Convert from Philips data/list format.')
        parser_p_dl.add_argument('data', help='.data file', type=Path)
        parser_p_dl.add_argument('list', help='.list file', type=Path)
        parser_p_dl.add_argument('spar', help='.SPAR file', type=Path)
        # for idx in range(5, 8):
        #     parser_p_dl.add_argument(f"-d{idx}", f"--dim{idx}", type=str, help=f"Specify dim {idx} loop counter.")
        #     parser_p_dl.add_argument(f"-t{idx}", f"--tag{idx}", type=str, help=f"Specify dim {idx} NIfTI MRS tag.")
        parser_p_dl = add_common_parameters(parser_p_dl)
        parser_p_dl.set_defaults(func=self.philips_dl)

        # Handle philips DICOM subcommand
        parser_philips_dicom = subparsers.add_parser('philips_dcm', help='Convert from Philips DICOM format.')
        parser_philips_dicom.add_argument('file', help='file or directory to convert', type=str)
        parser_philips_dicom.add_argument("-t", "--tag", type=str, help="Specify NIfTI MRS tag used for 5th "
                                                                        "dimension if multiple files are passed.")
        parser_philips_dicom = add_common_parameters(parser_philips_dicom)
        parser_philips_dicom.set_defaults(func=self.philips_dicom)

        # Handle GE subcommand
        parser_ge = subparsers.add_parser('ge', help='Convert from GE p-file format.')
        parser_ge.add_argument('file', help='file to convert', type=Path)
        parser_ge = add_common_parameters(parser_ge)
        parser_ge.set_defaults(func=self.ge)

        # # Handle ismrmrd subcommand
        # parser_ismrmrd = subparsers.add_parser('ismrmrd', help='Convert from ismrmrd format.')
        # parser_ismrmrd.add_argument('file', help='file to convert', type=str)
        # parser_ismrmrd.add_argument('a', type=str, help='placeholder')
        # parser_ismrmrd.set_defaults(func=self.ismrmrd)

        # Handle text subcommand
        parser_txt = subparsers.add_parser('text', help='Convert from plain text format.')
        parser_txt.add_argument('file', help='file to convert', type=str)
        parser_txt.add_argument("-i", "--imagingfreq", type=float,
                                help="Imaging (central) frequency in MHz", required=True)
        parser_txt.add_argument("-b", "--bandwidth", type=float,
                                help="Reciever bandwidth (sweepwidth) in Hz.", required=True)
        parser_txt.add_argument("-n", "--nucleus", type=str,
                                help="Nucleus string. e.g. 1H or 31P.", required=True)
        parser_txt.add_argument("-a", "--affine", type=str, help="NIfTI affine file", required=False, metavar='<file>')
        parser_txt = add_common_parameters(parser_txt)
        parser_txt.set_defaults(func=self.text)

        # jMRUI formats (.txt and .mrui)
        parser_jmrui = subparsers.add_parser('jmrui', help='Convert from jMRUI text or binary formats.')
        parser_jmrui.add_argument('file', help='file to convert', type=Path)
        parser_jmrui.add_argument("-a", "--affine", type=str, help="NIfTI affine file", required=False,
                                  metavar='<file>')
        parser_jmrui = add_common_parameters(parser_jmrui)
        parser_jmrui.set_defaults(func=self.jmrui)

        # LCModel RAW/H2O formats
        parser_raw = subparsers.add_parser('raw', help='Convert from LCModel RAW text format.')
        parser_raw.add_argument('file', help='file to convert', type=str)
        parser_raw.add_argument("-n", "--nucleus", type=str,
                                help="Nucleus string. e.g. 1H or 31P.", required=True)
        parser_raw.add_argument("-a", "--affine", type=str, help="NIfTI affine file", required=False, metavar='<file>')
        parser_raw = add_common_parameters(parser_raw)
        parser_raw.set_defaults(func=self.raw)

        # Bruker format
        parser_bruker = subparsers.add_parser('bruker', help='Convert from Bruker data format.')
        parser_bruker.add_argument('file', help='2dseq file to convert', type=str)
        parser_bruker.add_argument("-a", "--affine", type=str, help="NIfTI affine file", required=False, metavar='<file>')
        parser_bruker.add_argument('-q', '--query', action='append', default=[])
        parser_bruker = add_common_parameters(parser_bruker)
        parser_bruker.set_defaults(func=self.bruker)

        if len(sys.argv) == 1:
            parser.print_usage(sys.stderr)
            sys.exit(1)

        args = parser.parse_args()

        self.fileoutNames = []
        self.imageOut = []
        self.orientationInfoOut = []
        self.dwellTimes = []
        self.metaData = []

        self.outputDir = args.outdir

        if args.nifti1:
            nifti_mrs.NIfTI_MRS = nifti_mrs.get_mrs_class(nifti=1)
        else:
            nifti_mrs.NIfTI_MRS = nifti_mrs.get_mrs_class(nifti=2)

        args.func(args)

        if self.imageOut:
            self.validate_output()
            self.write_output(args.json)
        elif hasattr(args, 'view') and not args.view:
            print('No files to write.')

    def validate_output(self):
        """Run NIfTI MRS validation on output."""
        for f_out, nifti_mrs_img in zip(self.fileoutNames, self.imageOut):
            nifti_mrs_img.validate()

    def write_output(self, write_json=False):
        """Write any NIfTI MRS objects stored.
        If write_json is true also write meta-data as sidecar.
        """
        self.outputDir.mkdir(parents=True, exist_ok=True)

        for f_out, nifti_mrs_img in zip(self.fileoutNames, self.imageOut):
            out = self.outputDir / (f_out + '.nii.gz')
            nifti_mrs_img.save(out)

            if write_json:
                out_json = self.outputDir / (f_out + '.json')
                with open(out_json, 'w') as fp:
                    json.dump(json.loads(nifti_mrs_img.header.extensions[0].get_content()), fp, indent=4)

    # Start of the specific format handling functions.
    # Siemens twix (.dat) format
    def twix(self, args):
        """ Twix format handler."""
        from mapvbvd import mapVBVD
        from spec2nii.twixfunctions import process_twix, examineTwix

        # Call mapvbvd to load the twix file.
        twixObj = mapVBVD(args.file, quiet=args.quiet)

        if args.view:
            examineTwix(twixObj, op.basename(args.file), args.multiraid)
        else:
            if isinstance(twixObj, list):
                twixObj = twixObj[args.multiraid - 1]

            if not args.quiet:
                print(f"Converting twix file {args.file}.")
                print(f'Looking for evalinfo flag {args.evalinfo}.')

            overrides = {'dims': (args.dim5, args.dim6, args.dim7),
                         'tags': (args.tag5, args.tag6, args.tag7)}

            self.imageOut, self.fileoutNames = process_twix(twixObj,
                                                            args.fileout,
                                                            args.file.name,
                                                            args.evalinfo,
                                                            overrides,
                                                            args.quiet,
                                                            args.verbose)

    # (Siemens) DICOM (.ima) format
    def dicom(self, args):
        """Siemens DICOM format handler."""
        from spec2nii.dicomfunctions import multi_file_dicom
        path_in = Path(args.file)
        if path_in.is_dir():
            # Look for typical dicom file extensions
            files_in = sorted(path_in.glob('*.IMA')) + \
                sorted(path_in.glob('*.ima')) + \
                sorted(path_in.glob('*.dcm'))

            # If none found look for all files
            if len(files_in) == 0:
                files_in = sorted([x for x in path_in.iterdir() if x.is_file()])

            print(f'Found {len(files_in)} files.')
        else:
            print('Single file conversion.')
            files_in = [path_in]

        self.imageOut, self.fileoutNames = multi_file_dicom(files_in, args.fileout, args.tag, args.verbose)

    # (UIH) DICOM (.dcm) format
    def uih_dicom(self, args):
        """UIH DICOM format handler."""
        from spec2nii.uih import multi_file_dicom
        path_in = Path(args.file)
        if path_in.is_dir():
            # Look for typical dicom file extensions
            files_in = sorted(path_in.glob('*.dcm'))

            # If none found look for all files
            if len(files_in) == 0:
                files_in = sorted([x for x in path_in.iterdir() if x.is_file()])

            print(f'Found {len(files_in)} files.')
        else:
            print('Single file conversion.')
            files_in = [path_in]

        self.imageOut, self.fileoutNames = multi_file_dicom(files_in, args.fileout, args.tag, args.verbose)

    # Philips SDAP/SPAR handler
    def philips(self, args):
        # philips specific imports
        from spec2nii.philips import read_sdat_spar_pair

        self.imageOut = read_sdat_spar_pair(args.sdat, args.spar, args.tag)

        # name of output
        if args.fileout:
            mainStr = args.fileout
        elif args.sdat.stem == args.spar.stem:
            mainStr = args.sdat.stem
        else:
            mainStr = args.sdat.stem

        self.fileoutNames.append(mainStr)

    # Philips data/list (+SPAR) handler
    def philips_dl(self, args):
        # philips specific imports
        from spec2nii.philips_data_list import read_data_list_pair

        # TO DO
        # overrides = {'dims': (args.dim5, args.dim6, args.dim7),
        #              'tags': (args.tag5, args.tag6, args.tag7)}

        # self.imageOut, file_names = read_data_list_pair(args.data, args.list, args.spar, overrides)
        self.imageOut, file_names = read_data_list_pair(args.data, args.list, args.spar)

        # name of output
        if args.fileout:
            mainStr = args.fileout
        elif args.data.stem == args.list.stem:
            mainStr = args.data.stem
        else:
            mainStr = args.data.stem

        for fn in file_names:
            self.fileoutNames.append(mainStr + '_' + fn)

    # Philips DICOM handler
    def philips_dicom(self, args):
        """Philips DICOM format handler."""
        from warnings import warn
        from spec2nii.philips_dcm import multi_file_dicom
        warn('This Philips DICOM conversion routine is experimental and poorly tested.'
             ' Please get in contact with test data to help improve it.')
        path_in = Path(args.file)
        if path_in.is_dir():
            # Look for typical dicom file extensions
            files_in = sorted(path_in.glob('*.dcm'))

            # If none found look for all files
            if len(files_in) == 0:
                files_in = sorted([x for x in path_in.iterdir() if x.is_file()])

            print(f'Found {len(files_in)} files.')
        else:
            print('Single file conversion.')
            files_in = [path_in]

        self.imageOut, self.fileoutNames = multi_file_dicom(files_in, args.fileout, args.tag, args.verbose)

    # GE p file handler
    def ge(self, args):
        # ge specific imports
        from spec2nii.GE.ge_pfile import read_pfile
        self.imageOut, self.fileoutNames = read_pfile(args.file, args.fileout)

    def text(self, args):
        from spec2nii.other_formats import text
        # Simply pass through to the specific text function
        self.imageOut, self.fileoutNames = text(args)

    def jmrui(self, args):
        from spec2nii.jmrui import jmrui_format
        # Pass straight through to dedicated function
        self.imageOut, self.fileoutNames = jmrui_format(args)

    def raw(self, args):
        from spec2nii.other_formats import lcm_raw
        self.imageOut, self.fileoutNames = lcm_raw(args)

    # Bruker 2dseq files with FG_COMPLEX
    def bruker(self, args):
        from spec2nii.bruker import read_bruker
        self.imageOut, self.fileoutNames = read_bruker(args)


def main(*args):
    spec2nii(*args)
    return 0
