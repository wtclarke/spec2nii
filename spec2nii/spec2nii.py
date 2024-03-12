""" spec2nii - tool for conversion of various MRS data formats to NIFTI format.

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
from spec2nii import __version__ as spec2nii_ver
from numpy import isclose
# There are case specific imports below


class Spec2niiError(Exception):
    pass


class spec2nii:
    def __init__(self):
        cite_str = "Clarke WT, Bell TK, Emir UE, Mikkelsen M, Oeltzschner G, Shamaei A, Soher BJ, Wilson M. "\
                   "NIfTI-MRS: A standard data format for magnetic resonance spectroscopy. "\
                   "Magn Reson Med. 2022. doi: 10.1002/mrm.29418."
        parser = argparse.ArgumentParser(
            description='Convert raw spectroscopy data to NIfTI format.',
            epilog=f"If you use spec2nii please cite: {cite_str}")
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
            subparser.add_argument(
                '--anon',
                action='store_true',
                help="Create file without sensitive metadata. For greater control use spec2nii anon.")
            subparser.add_argument('--verbose', action='store_true')

            return subparser

        # Auto subcommand - heuristic ID of file type
        parser_auto = subparsers.add_parser('auto', help='Attempt automatic identification and conversion.')
        parser_auto.add_argument('file', help='file to convert', type=Path)
        parser_auto = add_common_parameters(parser_auto)
        parser_auto.set_defaults(func=self.auto)

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
        parser_rda = subparsers.add_parser('rda', help='Convert from Siemens spectroscopy .rda format.')
        parser_rda.add_argument('file', help='file to convert', type=Path)
        parser_rda = add_common_parameters(parser_rda)
        parser_rda.set_defaults(func=self.rda)

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
            "-t", "--tags",
            type=str,
            nargs='+',
            default=[None, None, None],
            help="Specify NIfTI MRS tags used for higher (5th-7th) dimensions. "
                 "Defaults to DIM_DYN if more than one spectrum is present. "
                 "Can be used to create singleton higher dimensions.")
        parser_philips.add_argument(
            "-s", "--shape",
            type=int,
            nargs='+',
            default=None,
            help="Specify the shape of higher (5th-7th) dimensions. Applies numpy style reshaping.")
        parser_philips.add_argument(
            "--special",
            type=str,
            default=None,
            help="Identify special case sequence. Options: 'hyper', 'hyper-ref'.")
        parser_philips = add_common_parameters(parser_philips)
        parser_philips.set_defaults(func=self.philips)

        # Handle philips data/list subcommand
        parser_p_dl = subparsers.add_parser('philips_dl', help='Convert from Philips data/list format.')
        parser_p_dl.add_argument('data', help='.data file', type=Path)
        parser_p_dl.add_argument('list', help='.list file', type=Path)
        parser_p_dl.add_argument('aux', help='Auxiliary file: .SPAR or DICOM file', type=Path)
        # for idx in range(5, 8):
        #     parser_p_dl.add_argument(f"-d{idx}", f"--dim{idx}", type=str, help=f"Specify dim {idx} loop counter.")
        #     parser_p_dl.add_argument(f"-t{idx}", f"--tag{idx}", type=str, help=f"Specify dim {idx} NIfTI MRS tag.")
        parser_p_dl.add_argument(
            "--special",
            type=str,
            default=None,
            help="Identify special case sequence. Options: 'hyper'.")
        parser_p_dl = add_common_parameters(parser_p_dl)
        parser_p_dl.set_defaults(func=self.philips_dl)

        # Handle philips DICOM subcommand
        parser_philips_dicom = subparsers.add_parser('philips_dcm', help='Convert from Philips DICOM format.')
        parser_philips_dicom.add_argument('file', help='file or directory to convert', type=str)
        parser_philips_dicom.add_argument(
            "-t",
            "--tag",
            type=str,
            default='DIM_DYN',
            help="Specify NIfTI MRS tag used for 5th dimension if multiple files are passed. "
                 "Defaults to DIM_DYN.")
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
        parser_raw.add_argument("-i", "--imagingfreq", type=float,
                                help="Optional. Imaging (central) frequency in MHz", required=False)
        parser_raw.add_argument("-b", "--bandwidth", type=float,
                                help="Optional. Receiver bandwidth (spectral width) in Hz.", required=False)
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
                                 override_dwelltime=None,
                                 anon=False)

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
        parser_insert.add_argument("--dwelltime", type=float,
                                   help="Specify a new dwelltime (1/bandwidth, pixdim[4]) value in seconds.")
        parser_insert.add_argument("-f", "--fileout", type=str,
                                   help="Output file base name (default = input file name)")
        parser_insert.add_argument("-o", "--outdir", type=Path,
                                   help="Output location (default = .)", default='.')
        parser_insert.set_defaults(func=self.insert,
                                   nifti1=False,
                                   json=False,
                                   override_nucleus=None,
                                   override_frequency=None,
                                   override_dwelltime=None,
                                   verbose=False,
                                   anon=False)

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

        args.func(args)

        if self.imageOut:
            self.implement_overrides(args)

            self.insert_spectralwidth()

            if args.anon:
                from spec2nii.anonymise import anon_nifti_mrs
                for idx, nifti_mrs_img in enumerate(self.imageOut):
                    self.imageOut[idx] = anon_nifti_mrs(nifti_mrs_img, verbose=args.verbose)

            self.validate_output()
            self.write_output(args.json, args.nifti1)
            self.validate_write(args.verbose)
            if args.verbose:
                print(f'Please cite {cite_str}.')
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

    def insert_spectralwidth(self):
        """Ensure that the correct spectral width is inserted into the header extension"""
        for nifti_mrs_img in self.imageOut:
            if 'SpectralWidth' in nifti_mrs_img.hdr_ext\
                    and not isclose(
                        nifti_mrs_img.hdr_ext['SpectralWidth'],
                        1 / nifti_mrs_img.dwelltime,
                        atol=1E-2):
                nifti_mrs_img.remove_hdr_field('SpectralWidth')
                nifti_mrs_img.add_hdr_field('SpectralWidth', 1 / nifti_mrs_img.dwelltime)
            else:
                nifti_mrs_img.add_hdr_field('SpectralWidth', 1 / nifti_mrs_img.dwelltime)

    def validate_output(self):
        """Run NIfTI MRS validation on output."""
        import nifti_mrs.validator as validate
        # Currently this repeats the validation before the save.
        # But useful here to do exception handling
        for f_out, nifti_mrs_img in zip(self.fileoutNames, self.imageOut):
            try:
                validate.validate_nifti_mrs(nifti_mrs_img)
            except (
                    validate.headerExtensionError,
                    validate.niftiDataError,
                    validate.niftiHeaderError) as exc:

                raise Spec2niiError(f'Generated file {f_out} failed validation.') from exc

    def write_output(self, write_json=False, nifti1=False):
        """Write any NIfTI MRS objects stored.
        If write_json is true also write meta-data as sidecar.
        """
        self.outputDir.mkdir(parents=True, exist_ok=True)

        for f_out, nifti_mrs_img in zip(self.fileoutNames, self.imageOut):
            out = self.outputDir / (f_out + '.nii.gz')

            if nifti1:
                # If nifti1 is requested just remake the file.
                # This is more than a bit hacky, but avoids passing the option to every end-point.
                from nifti_mrs.create_nmrs import gen_nifti_mrs_hdr_ext
                gen_nifti_mrs_hdr_ext(
                    nifti_mrs_img[:],
                    nifti_mrs_img.dwelltime,
                    nifti_mrs_img.hdr_ext,
                    nifti_mrs_img.getAffine('voxel', 'world'),
                    nifti_version=1)\
                    .save(out)
            else:
                nifti_mrs_img.save(out)

            if write_json:
                out_json = self.outputDir / (f_out + '.json')
                with open(out_json, 'w') as fp:
                    json.dump(json.loads(nifti_mrs_img.header.extensions[0].get_content()), fp, indent=4)

    class NIfTIMRSWriteError(IOError):
        pass

    def validate_write(self, verbose):
        for f_out in self.fileoutNames:
            out = self.outputDir / (f_out + '.nii.gz')
            if out.exists() and verbose:
                print(f'Output {out.name} written to {out.parent}')
            elif out.exists():
                pass
            else:
                raise self.NIfTIMRSWriteError(f'Output {out.name} in {out.parent} not found!')

    # Attempt automatic file type identification
    def auto(self, args):
        """Attempt automatic file type identification and conversion

        Currently handles: Twix, RDA, SPAR/SDAT, GE pfile, DICOM
        """
        from warnings import warn
        warn('Automatic conversion is an experimental feature. Please verify the output carefully.')

        # Twix format - ID by extension
        if args.file.suffix.lower() == '.dat':
            print("Attempting conversion as Siemens Twix (.dat), evalinfo = 'image'")
            setattr(args, 'evalinfo', 'image')
            setattr(args, 'quiet', False)
            setattr(args, 'view', False)
            for idx in range(5, 8):
                setattr(args, f'dim{idx}', None)
                setattr(args, f'tag{idx}', None)
            setattr(args, 'remove_os', False)
            self.twix(args)

        # Siemens RDA format - ID by extension
        elif args.file.suffix.lower() == '.rda':
            print('Attempting conversion as Siemens RDA (.rda)')
            self.rda(args)

        # Philips SPAR/SDAT format - ID by extension(s)
        elif args.file.suffix.lower() in ('.spar', '.sdat'):
            print('Attempting conversion as Philips SPAR/SDAT pair')
            if args.file.with_suffix('.SPAR').is_file()\
                    and args.file.with_suffix('.SDAT').is_file():
                setattr(args, 'spar', args.file.with_suffix('.SPAR'))
                setattr(args, 'sdat', args.file.with_suffix('.SDAT'))
            else:
                raise Spec2niiError(
                    'Unable to find both files in SPAR/SDAT pair.'
                    'Please use manual selection of subcomands. '
                    'e.g. spec2nii philips ...')
            setattr(args, 'tags', ["DIM_DYN", None, None])
            setattr(args, 'shape', None)
            setattr(args, 'special', None)
            self.philips(args)

        # Philips DATA/LIST - Reject because of unknown aux file format.
        elif args.file.suffix.lower() in ('.data', '.list'):
            raise Spec2niiError(
                'Automatic conversion not setup for data/list conversion. '
                'Please use manual selection of subcomand. '
                'e.g. spec2nii philips_dl ...')

        # GE pfile (.7) format - ID by extension
        elif args.file.suffix.lower() == '.7':
            self.ge(args)

        # If no matches assume DICOM - ID by loading file
        else:
            import pydicom as pdcm
            try:
                if args.file.is_dir():
                    files_in = \
                        sorted(args.file.rglob('*.IMA')) + \
                        sorted(args.file.rglob('*.ima')) + \
                        sorted(args.file.rglob('*.dcm'))
                    file = pdcm.read_file(files_in[0])
                else:
                    file = pdcm.read_file(args.file)

                manufacturer = file.Manufacturer
                setattr(args, 'tag', None)
                if manufacturer.lower() == 'siemens':
                    print('Attempting conversion as Siemens DICOM.')
                    setattr(args, 'voi', False)
                    self.dicom(args)
                elif manufacturer.lower() == 'uih':
                    print('Attempting conversion as UIH DICOM.')
                    self.uih_dicom(args)
                elif manufacturer.lower() == 'philips medical systems':
                    print('Attempting conversion as Philips DICOM.')
                    self.philips_dicom(args)
                else:
                    raise Spec2niiError(f'Unknown DICOM manufacturer {manufacturer}.')
            # No successful ID as DICOM - fail at automatic load.
            except pdcm.errors.InvalidDicomError:
                raise Spec2niiError(
                    'Unable to automatically identify file type. '
                    'Please use manual selection of subcomands. '
                    'e.g. spec2nii twix ..., spec2nii philips ...')

    # Start of the specific format handling functions.
    # Siemens twix (.dat) format
    def twix(self, args):
        """ Twix format handler."""
        from mapvbvd import mapVBVD
        from spec2nii.Siemens.twixfunctions import process_twix, examineTwix

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
        from spec2nii.Siemens.dicomfunctions import multi_file_dicom
        path_in = Path(args.file)
        if path_in.is_dir():
            # Look for typical dicom file extensions
            files_in = sorted(path_in.rglob('*.IMA')) + \
                sorted(path_in.rglob('*.ima')) + \
                sorted(path_in.rglob('*.dcm'))

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
                mask_data = np.ones((1, 1, 1), int)
                save(
                    Nifti2Image(
                        mask_data,
                        img.hdr_ext['VOI'],
                        dtype=mask_data.dtype),
                    out)

    # Siemens RDA (.rda) format
    def rda(self, args):
        """Siemens RDA format handler."""
        from spec2nii.Siemens.rda import convert_rda
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
        from spec2nii.Philips.philips import read_sdat_spar_pair

        self.imageOut, self.fileoutNames = read_sdat_spar_pair(
            args.sdat,
            args.spar,
            args.shape,
            args.tags,
            args.fileout,
            args.special)

    # Philips data/list (+SPAR) handler
    def philips_dl(self, args):
        # philips specific imports
        from spec2nii.Philips.philips_data_list import read_data_list_pair

        # TO DO
        # overrides = {'dims': (args.dim5, args.dim6, args.dim7),
        #              'tags': (args.tag5, args.tag6, args.tag7)}

        # self.imageOut, file_names = read_data_list_pair(args.data, args.list, args.spar, overrides)
        self.imageOut, file_names = read_data_list_pair(args.data, args.list, args.aux, args.special)

        # name of output
        if args.fileout:
            mainStr = args.fileout
        else:
            mainStr = args.data.stem

        for fn in file_names:
            self.fileoutNames.append(mainStr + '_' + fn)

    # Philips DICOM handler
    def philips_dicom(self, args):
        """Philips DICOM format handler."""
        from warnings import warn
        from spec2nii.Philips.philips_dcm import multi_file_dicom
        warn('This Philips DICOM conversion routine has limited testing outside vendor PRESS and MEGA sequences.'
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
        from spec2nii.anonymise import anon_file
        self.imageOut, self.fileoutNames = anon_file(args)

    # Dump function
    @staticmethod
    def dump(args):
        from spec2nii.other import dump_headers
        dump_headers(args)

    # Extract function
    @staticmethod
    def extract(args):
        from spec2nii.other import extract_hdr_ext
        extract_hdr_ext(args)

    # Insert function
    def insert(self, args):
        from spec2nii.other import insert_hdr_ext
        self.imageOut, self.fileoutNames = insert_hdr_ext(args)


def main(*args):
    spec2nii(*args)
    return 0
