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

from nibabel.nifti2 import Nifti2Image
from spec2nii import nifti_mrs
from spec2nii import __version__ as spec2nii_ver
# There are case specific imports below


class spec2nii:
    def __init__(self):
        parser = argparse.ArgumentParser(description='Convert raw spectroscopy data to NIfTI format.')
        parser.add_argument('-v', '--version', action='version', version=spec2nii_ver)

        subparsers = parser.add_subparsers(title='spec2nii subcommands')

        def add_common_parameters(subparser):
            # Add options that are common to all subcommands
            subparser.add_argument('-j', '--json', help='Create json sidecar.', action='store_true')
            subparser.add_argument("-f", "--fileout", type=str,
                                   help="Output file base name (default = input file name)")
            subparser.add_argument("-o", "--outdir", type=Path,
                                   help="Output location (default = .)", default='.')
            subparser.add_argument('--nifti1', action='store_true')
            subparser.add_argument("--override_nucleus", type=str, nargs='+',
                                   help="Override ResonantNucleus field with input(s). E.g. '2H'.")
            subparser.add_argument("--override_frequency", type=float, nargs='+',
                                   help="Override SpectrometerFrequency field with input(s). Input in MHz.")
            subparser.add_argument("--override_dwelltime", type=float,
                                   help="Override dwell time field with input. Input in seconds.")
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
        parser_twix.add_argument('--remove_os', action='store_true', help='Remove time-domain oversampling.')
        parser_twix = add_common_parameters(parser_twix)
        parser_twix.set_defaults(func=self.twix)

        # Handle dicom subcommand
        parser_dicom = subparsers.add_parser('dicom', help='Convert from Siemens DICOM format.')
        parser_dicom.add_argument('file', help='file or directory to convert', type=str)
        parser_dicom.add_argument("-t", "--tag", type=str, help="Specify NIfTI MRS tag used for 5th "
                                                                "dimension if multiple files are passed.")
        parser_dicom.add_argument('--voi', action='store_true', help='Output VOI as single voxel NIfTI mask.')
        parser_dicom = add_common_parameters(parser_dicom)
        parser_dicom.set_defaults(func=self.dicom)

        # Handle rda subcommand
        parser_dicom = subparsers.add_parser('rda', help='Convert from Siemens spectroscopy .rda format.')
        parser_dicom.add_argument('file', help='file to convert', type=Path)
        parser_dicom = add_common_parameters(parser_dicom)
        parser_dicom.set_defaults(func=self.rda)

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
        parser_philips.add_argument(
            "-t", "--tag",
            type=str,
            default="DIM_DYN",
            help="Specify NIfTI MRS tag used for 5th dimension if multiple transients are present.")
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
                                help="Receiver bandwidth (sweepwidth) in Hz.", required=True)
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
        parser_bruker.add_argument('-q', '--query', action='append', default=[])
        parser_bruker.add_argument('-m', '--mode', type=str, default='2DSEQ', choices=['2DSEQ', 'FID'])
        parser_bruker.add_argument('-d', '--dump_headers',
                                   help='Dump bruker header files into json header extension',
                                   action='store_true')
        parser_bruker = add_common_parameters(parser_bruker)
        parser_bruker.set_defaults(func=self.bruker)

        # Varian format

        parser_varian = subparsers.add_parser('varian', help='Convert from the Varian/OpenVnmrJ data format')
        parser_varian.add_argument(
            'file', help='Path to the varian .fid directory, containing procpar and fid files', type=str)
        parser_varian.add_argument("-t", "--tag6", type=str, default='DIM_DYN',
                                   help="Specify 6th dimension NIfTI MRS tag (default DIM_DYN).")
        parser_varian.add_argument('-d', '--dump_headers',
                                   help='Dump varian header information into json header extension',
                                   action='store_true')
        parser_varian = add_common_parameters(parser_varian)
        parser_varian.set_defaults(func=self.varian)

        # Additional functions - anonymise and dump
        parser_anon = subparsers.add_parser('anon', help='Anonymise existing NIfTI-MRS file.')
        parser_anon.add_argument('file', help='file to anonymise', type=Path)
        parser_anon.add_argument("-r", "--remove", type=str, metavar='KEY',
                                 action='append',
                                 help="Explicitly remove key. Argument may be repeated")
        parser_anon.add_argument("-f", "--fileout", type=str,
                                 help="Output file base name (default = input file name)")
        parser_anon.add_argument("-o", "--outdir", type=Path,
                                 help="Output location (default = .)", default='.')
        parser_anon.add_argument('-v', '--verbose', action='store_true')
        parser_anon.set_defaults(func=self.anon,
                                 json=False,
                                 nifti1=False,
                                 override_nucleus=None,
                                 override_frequency=None,
                                 override_dwelltime=None)

        parser_dump = subparsers.add_parser('dump', help='Dump contents of headers from existing NIfTI-MRS file.')
        parser_dump.add_argument('file', help='NIfTI-MRS file', type=Path)
        parser_dump.set_defaults(func=self.dump,
                                 outdir=None,
                                 nifti1=False)

        parser_extract = subparsers.add_parser('extract', help='Generate json sidecar from existing NIfTI-MRS file.')
        parser_extract.add_argument('file', help='NIFTI-MRS file', type=Path)
        parser_extract.add_argument("-o", "--outdir", type=Path,
                                    help="Output location (default: directory of input file)", default=None)
        parser_extract.set_defaults(func=self.extract,
                                    nifti1=False)

        parser_insert = subparsers.add_parser('insert', help='Insert json formatted file into existing NIfTI-MRS file.')
        parser_insert.add_argument('file', help='NIFTI-MRS file', type=Path)
        parser_insert.add_argument('json_file', help='JSON file to insert', type=Path)
        parser_insert.add_argument("-f", "--fileout", type=str,
                                   help="Output file base name (default = input file name)")
        parser_insert.add_argument("-o", "--outdir", type=Path,
                                   help="Output location (default = .)", default='.')
        parser_insert.set_defaults(func=self.insert,
                                   nifti1=False,
                                   json=False,
                                   override_nucleus=None,
                                   override_frequency=None,
                                   override_dwelltime=None)

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
            self.implement_overrides(args)
            self.validate_output()
            self.write_output(args.json)
        elif hasattr(args, 'view') and not args.view:
            print('No files to write.')

    def implement_overrides(self, args):
        """Implement any command line overrides for essential parameters."""
        for nifti_mrs_img in self.imageOut:
            if args.override_dwelltime:
                nifti_mrs_img.set_dwell_time(args.override_dwelltime)

            if args.override_nucleus or args.override_frequency:
                from nibabel.nifti1 import Nifti1Extension
                hdr_ext_codes = nifti_mrs_img.header.extensions.get_codes()
                index = hdr_ext_codes.index(44)
                original = json.loads(nifti_mrs_img.header.extensions[index].get_content())

                if args.override_nucleus:
                    original['ResonantNucleus'] = args.override_nucleus
                if args.override_frequency:
                    original['SpectrometerFrequency'] = args.override_frequency
                json_s = json.dumps(original)
                new_ext = Nifti1Extension(44, json_s.encode('UTF-8'))
                nifti_mrs_img.header.extensions.clear()
                nifti_mrs_img.header.extensions.append(new_ext)

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
                                                            args.verbose,
                                                            remove_os=args.remove_os)

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

        # VOI mask - TO DO extend to other formats
        if args.voi:
            import numpy as np
            from nibabel import save
            for img, f_out in zip(self.imageOut, self.fileoutNames):
                out = self.outputDir / (f_out + '_voi.nii.gz')
                save(Nifti2Image(np.ones((1, 1, 1), int), img.hdr_ext['VOI']), out)

    # Siemens RDA (.rda) format
    def rda(self, args):
        """Siemens RDA format handler."""
        from spec2nii.rda import convert_rda
        self.imageOut, self.fileoutNames = convert_rda(args.file, args.fileout, args.verbose)

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

    # Varian parser
    def varian(self, args):
        from spec2nii.varian_importer import read_varian
        self.imageOut, self.fileoutNames = read_varian(args)

    # Anonymise function
    def anon(self, args):
        from spec2nii.anonymise import anon_nifti_mrs
        self.imageOut, self.fileoutNames = anon_nifti_mrs(args)

    # Dump function
    def dump(self, args):
        from spec2nii.other import dump_headers
        dump_headers(args)

    # Extract function
    def extract(self, args):
        from spec2nii.other import extract_hdr_ext
        extract_hdr_ext(args)

    # Insert function
    def insert(self, args):
        from spec2nii.other import insert_hdr_ext
        self.imageOut, self.fileoutNames = insert_hdr_ext(args)


def main(*args):
    spec2nii(*args)
    return 0
