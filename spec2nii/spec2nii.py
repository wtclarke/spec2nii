""" spec2nii - tool for conversion of various MRS data formats to NIFTI format.

spec2nii converts the following formats to NIFTI files.
Supporting SVS:
    Siemens "Twix" .dat format
    Siemens DICOM
    Philips SPAR/SDAT files
    GE p-files
    LCModel RAW
    jMRUI text
    Plain text

Supporting CSI/MRSI:
    Siemens DICOM

This module contains the main class to be called as a script (through the main function).

Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford
"""

import argparse
import sys
import numpy as np
import os.path as op
from pathlib import Path
from spec2nii.nifti_orientation import NIFTIOrient
import json
# There are case specific imports below


class spec2nii:
    def __init__(self):
        parser = argparse.ArgumentParser(description='Convert raw spectroscopy data to NIfTI format.')

        subparsers = parser.add_subparsers(title='subcommands',
                                           description='File types supported')

        # Handle twix subcommand
        parser_twix = subparsers.add_parser('twix', help='Convert from Siemens .dat twix format.')
        parser_twix.add_argument('file', help='file to convert', type=Path)
        parser_twix.add_argument('-j', '--json', help='Create json sidecar.', action='store_true')
        group = parser_twix.add_mutually_exclusive_group(required=True)
        group.add_argument("-v", "--view", help="View contents of twix file, no files converted", action='store_true')
        group.add_argument('-e', '--evalinfo', type=str, help='evalInfo flag to convert')
        parser_twix.add_argument("-f", "--fileout", type=str, help="Output file base name (default = input file name)")
        parser_twix.add_argument("-o", "--outdir", type=Path, help="Output location (default = .)", default='.')
        parser_twix.add_argument("-m", "--multiraid", type=int, help="Select multiraid file to load"
                                 " (default = 2 i.e. 2nd file)", default=2)
        parser_twix.add_argument("-q", "--quiet", help="Suppress text output", action='store_true')
        for idx in range(5, 8):
            parser_twix.add_argument(f"-d{idx}", f"--dim{idx}", type=str, help=f"Specify dim {idx} loop counter.")
            parser_twix.add_argument(f"-t{idx}", f"--tag{idx}", type=str, help=f"Specify dim {idx} NIfTI MRS tag.")
        parser_twix.set_defaults(func=self.twix)

        # Handle dicom subcommand
        parser_dicom = subparsers.add_parser('dicom', help='Convert from DICOM format.')
        parser_dicom.add_argument('file', help='file or directory to convert', type=str)
        parser_dicom.add_argument("-f", "--fileout", type=str,
                                  help="Output file base name (default = DCM SeriesDescription tag)")
        parser_dicom.add_argument("-o", "--outdir", type=str, help="Output location (default = .)", default='.')
        parser_dicom.add_argument('-j', '--json', help='Create json sidecar.', action='store_true')
        parser_dicom.set_defaults(func=self.dicom)

        # Handle philips subcommand
        parser_philips = subparsers.add_parser('philips', help='Convert from Philips spar/sdat format.')
        parser_philips.add_argument('sdat', help='SDAT file', type=str)
        parser_philips.add_argument('spar', help='SPAR file', type=str)
        parser_philips.add_argument("-f", "--fileout", type=str,
                                    help="Output file base name (default = SDAT/SPAR name)")
        parser_philips.add_argument("-o", "--outdir", type=str, help="Output location (default = .)", default='.')
        parser_philips.add_argument('-j', '--json', help='Create json sidecar.', action='store_true')
        parser_philips.set_defaults(func=self.philips)

        # Handle GE subcommand
        parser_ge = subparsers.add_parser('ge', help='Convert from GE p-file format.')
        parser_ge.add_argument('file', help='file to convert', type=str)
        parser_ge.add_argument("-f", "--fileout", type=str, help="Output file name (default = input name)")
        parser_ge.add_argument("-o", "--outdir", type=str, help="Output location (default = .)", default='.')
        parser_ge.add_argument('-j', '--json', help='Create json sidecar.', action='store_true')
        parser_ge.set_defaults(func=self.ge)

        # # Handle ismrmrd subcommand
        # parser_ismrmrd = subparsers.add_parser('ismrmrd', help='Convert from ismrmrd format.')
        # parser_ismrmrd.add_argument('file', help='file to convert', type=str)
        # parser_ismrmrd.add_argument('a', type=str, help='placeholder')
        # parser_ismrmrd.set_defaults(func=self.ismrmrd)

        # Handle text subcommand
        parser_txt = subparsers.add_parser('text', help='Convert from plain text format.')
        parser_txt.add_argument('file', help='file to convert', type=str)
        parser_txt.add_argument('-j', '--json', help='Create json sidecar.', action='store_true')
        parser_txt.add_argument("-i", "--imagingfreq", type=float,
                                help="Imaging (central) frequency in MHz", required=True)
        parser_txt.add_argument("-b", "--bandwidth", type=float,
                                help="Reciever bandwidth (sweepwidth) in Hz.", required=True)
        parser_txt.add_argument("-a", "--affine", type=str, help="NIfTI affine file", required=False, metavar='<file>')
        parser_txt.add_argument("-f", "--fileout", type=str, help="Output file base name (default = input file name)")
        parser_txt.add_argument("-o", "--outdir", type=str, help="Output location (default = .)", default='.')
        parser_txt.set_defaults(func=self.text)

        parser_jmrui = subparsers.add_parser('jmrui', help='Convert from jMRUI text format.')
        parser_jmrui.add_argument('file', help='file to convert', type=str)
        parser_jmrui.add_argument('-j', '--json', help='Create json sidecar.', action='store_true')
        parser_jmrui.add_argument("-a", "--affine", type=str, help="NIfTI affine file", required=False,
                                  metavar='<file>')
        parser_jmrui.add_argument("-f", "--fileout", type=str, help="Output file base name (default = input file name)")
        parser_jmrui.add_argument("-o", "--outdir", type=str, help="Output location (default = .)", default='.')
        parser_jmrui.set_defaults(func=self.jmrui)

        parser_raw = subparsers.add_parser('raw', help='Convert from LCModel RAW text format.')
        parser_raw.add_argument('file', help='file to convert', type=str)
        parser_raw.add_argument('-j', '--json', help='Create json sidecar.', action='store_true')
        parser_raw.add_argument("-a", "--affine", type=str, help="NIfTI affine file", required=False, metavar='<file>')
        parser_raw.add_argument("-f", "--fileout", type=str, help="Output file base name (default = input file name)")
        parser_raw.add_argument("-o", "--outdir", type=str, help="Output location (default = .)", default='.')
        parser_raw.set_defaults(func=self.raw)

        parser.add_argument('--verbose', action='store_true')

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
        path_in = Path(args.file)
        if path_in.is_dir():
            # Look for typical dicom file extensions
            files_in = sorted(path_in.glob('*.IMA')) + \
                       sorted(path_in.glob('*.ima')) + \
                       sorted(path_in.glob('*.dcm'))

            # If none found look for all files
            if len(files_in) ==0:
                files_in = sorted([x for x in path_in.iterdir() if x.is_file()])

            print(f'Found {len(files_in)} files.')
        else:
            print('Single file conversion.')
            files_in = [path_in]
        
        # DICOM specific imports
        import nibabel.nicom.dicomwrappers
        from spec2nii.dicomfunctions import svs_or_CSI,process_siemens_svs,process_siemens_csi

        for idx,fn in enumerate(files_in):
            print(f'Converting dicom file {fn}')
            
            img = nibabel.nicom.dicomwrappers.wrapper_from_file(fn)

            mrs_type = svs_or_CSI(img)

            if mrs_type == 'SVS':
                specDataCmplx,currNiftiOrientation,dwelltime,meta = process_siemens_svs(img,args)
            
                if args.fileout:
                    mainStr = args.fileout
                else: 
                    mainStr = img.dcm_data.SeriesDescription
                self.fileoutNames.append(f'{mainStr}_{idx :03d}')

                newshape = (1,1,1)+specDataCmplx.shape
                specDataCmplx = specDataCmplx.reshape(newshape)

                self.imageOut.append(specDataCmplx)
                self.orientationInfoOut.append(currNiftiOrientation)
                self.dwellTimes.append(dwelltime)
                if args.json:
                    self.metaData.append(meta)

            elif mrs_type == 'CSI':
                specDataCmplx,currNiftiOrientation,dwelltime,meta = process_siemens_csi(img,args)
            
                if args.fileout:
                    mainStr = args.fileout
                else: 
                    mainStr = img.dcm_data.SeriesDescription
                self.fileoutNames.append(f'{mainStr}_{idx :03d}')

                self.imageOut.append(specDataCmplx)
                self.orientationInfoOut.append(currNiftiOrientation)
                self.dwellTimes.append(dwelltime)
                if args.json:
                    self.metaData.append(meta)


    def philips(self,args):
        # philips specific imports
        from spec2nii.philips import read_sdat_spar_pair

        data, orientation, dwelltime, meta = read_sdat_spar_pair(args.sdat,args.spar)

        # name of output
        if args.fileout:
            mainStr = args.fileout
        elif op.splitext(op.basename(args.sdat))[0]==op.splitext(op.basename(args.spar))[0]:
            mainStr = op.splitext(op.basename(args.sdat))[0]
        else:
            mainStr = op.splitext(op.basename(args.sdat))[0]

        newshape = (1,1,1)+data.shape
        data = data.reshape(newshape)

        # Place in data output format
        self.imageOut.append(data)
        self.orientationInfoOut.append(orientation)
        self.dwellTimes.append(dwelltime)
        self.metaData.append(meta)
        self.fileoutNames.append(mainStr)

    def ge(self,args):
        # philips specific imports
        from spec2nii.GE import read_p_file

        data,ref_data, orientation, dwelltime, meta = read_p_file(args.file)

        # name of output
        if args.fileout:
            baseStr = args.fileout
        else:
            baseStr = op.splitext(op.basename(args.file))[0]

        # Place in data output format
        for idx,d in enumerate(data):
            d = d.T
            newshape = (1,1,1)+d.shape
            d = d.reshape(newshape)
            self.imageOut.append(d)
            self.orientationInfoOut.append(orientation)
            self.dwellTimes.append(dwelltime)
            self.metaData.append(meta)
            self.fileoutNames.append(baseStr+f'_frame{idx:03.0f}')
        if ref_data is not None:
            for idx,d in enumerate(ref_data):
                d = d.T
                newshape = (1,1,1)+d.shape
                d = d.reshape(newshape)
                self.imageOut.append(d)
                self.orientationInfoOut.append(orientation)
                self.dwellTimes.append(dwelltime)
                self.metaData.append(meta)
                self.fileoutNames.append(baseStr+f'_ref_frame{idx:03.0f}')

    def text(self,args):
        # Read text from file
        data = np.loadtxt(args.file)
        data = data[:,0] + 1j*data[:,1] 
        
        newshape = (1,1,1)+data.shape
        data = data.reshape(newshape)

        # Interpret required arguments (frequency and bandwidth)
        imagingfreq = args.imagingfreq
        dwelltime = 1.0/args.bandwidth
        metadict = {'ImagingFrequency':imagingfreq,'Dwelltime':dwelltime}
        
        # Read optional affine file
        if args.affine:
            affine = np.loadtxt(args.affine)
        else:
            affine = np.eye(4)
        # qb,qc,qd,qx,qy,qz,dx,dy,dz,qfac = nifti_mat44_to_quatern(affine)        
        currNiftiOrientation = NIFTIOrient(affine)

        # File names
        if args.fileout:
            mainStr = args.fileout
        else:
            base=op.basename(args.file)
            mainStr = op.splitext(base)[0]

        # Place in data output format
        self.imageOut.append(data)
        self.orientationInfoOut.append(currNiftiOrientation)
        self.dwellTimes.append(dwelltime)
        self.metaData.append(metadict)
        self.fileoutNames.append(mainStr)

    def jmrui(self,args):
        from fsl_mrs.utils.mrs_io import jmrui_io
        # Read data from file
        data,header = jmrui_io.readjMRUItxt(args.file)

        newshape = (1,1,1)+data.shape
        data = data.reshape(newshape)

        # meta
        dwelltime = header['dwelltime']
        metadict = {'ImagingFrequency':header['centralFrequency'],'Dwelltime':header['dwelltime']}
        
        # Read optional affine file
        if args.affine:
            affine = np.loadtxt(args.affine)
        else:
            affine = np.eye(4)
        # qb,qc,qd,qx,qy,qz,dx,dy,dz,qfac = nifti_mat44_to_quatern(affine)        
        currNiftiOrientation = NIFTIOrient(affine)

        # File names
        if args.fileout:
            mainStr = args.fileout
        else:
            base=op.basename(args.file)
            mainStr = op.splitext(base)[0]

        # Place in data output format
        self.imageOut.append(data)
        self.orientationInfoOut.append(currNiftiOrientation)
        self.dwellTimes.append(dwelltime)
        self.metaData.append(metadict)
        self.fileoutNames.append(mainStr)

    def raw(self,args):
        from fsl_mrs.utils.mrs_io import lcm_io
        # Read data from file
        data,header = lcm_io.readLCModelRaw(args.file,conjugate=True)
        
        newshape = (1,1,1)+data.shape
        data = data.reshape(newshape)

        # meta
        dwelltime = header['dwelltime']
        metadict = {'ImagingFrequency':header['centralFrequency'],'Dwelltime':header['dwelltime']}
        
        # Read optional affine file
        if args.affine:
            affine = np.loadtxt(args.affine)
        else:
            affine = np.eye(4)
        # qb,qc,qd,qx,qy,qz,dx,dy,dz,qfac = nifti_mat44_to_quatern(affine)        
        currNiftiOrientation = NIFTIOrient(affine)

        # File names
        if args.fileout:
            mainStr = args.fileout
        else:
            base=op.basename(args.file)
            mainStr = op.splitext(base)[0]

        # Place in data output format
        self.imageOut.append(data)
        self.orientationInfoOut.append(currNiftiOrientation)
        self.dwellTimes.append(dwelltime)
        self.metaData.append(metadict)
        self.fileoutNames.append(mainStr)

def main(*args):
    spec2nii(*args)
    return 0

# if __name__== "__main__":
#     spec2nii()