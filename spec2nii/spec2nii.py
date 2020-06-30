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
    

This module contains the main class to be called as a script (through the main function).

Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford 
"""

import argparse
import sys
import numpy as np
from spec2nii.writeNii import writeNii
from spec2nii.writeJSON import writeJSON
import os.path as op
from pathlib import Path
from spec2nii.dcm2niiOrientation.orientationFuncs import nifti_dicom2mat
from spec2nii.nifti_orientation import NIFTIOrient
# There are case specific imports below


class spec2nii:
    def __init__(self):
        parser = argparse.ArgumentParser(description='Convert raw spectroscopy data to NIfTI format.')      

        subparsers = parser.add_subparsers(title='subcommands',
                                            description='File types supported')

        # Handle twix subcommand
        parser_twix = subparsers.add_parser('twix', help='Convert from Siemens .dat twix format.')
        parser_twix.add_argument('file',help='file to convert',type=str)
        parser_twix.add_argument('-j','--json',help='Create json sidecar.',action='store_true')
        group = parser_twix.add_mutually_exclusive_group(required=True)
        group.add_argument("-v", "--view", help="View contents of twix file, no files converted",action='store_true')
        group.add_argument('-e','--evalinfo', type=str, help='evalInfo flag to convert')
        parser_twix.add_argument("-f", "--fileout", type=str,help="Output file base name (default = input file name)")
        parser_twix.add_argument("-o", "--outdir", type=str,help="Output location (default = .)",default='.')
        parser_twix.add_argument("-m", "--multiraid", type=int,help="Select multiraid file to load (default = 2 i.e. 2nd file)",default=2)      
        parser_twix.add_argument("-q", "--quiet",help="Suppress text output",action='store_true')       
        parser_twix.set_defaults(func=self.twix)

        # Handle dicom subcommand
        parser_dicom = subparsers.add_parser('dicom', help='Convert from DICOM format.')
        parser_dicom.add_argument('file',help='file or directory to convert',type=str)
        parser_dicom.add_argument("-f", "--fileout", type=str,help="Output file base name (default = DCM SeriesDescription tag)")
        parser_dicom.add_argument("-o", "--outdir", type=str,help="Output location (default = .)",default='.')
        parser_dicom.add_argument('-j','--json',help='Create json sidecar.',action='store_true')
        parser_dicom.set_defaults(func=self.dicom)

        # Handle philips subcommand
        parser_philips = subparsers.add_parser('philips', help='Convert from Philips spar/sdat format.')
        parser_philips.add_argument('sdat',help='SDAT file',type=str)
        parser_philips.add_argument('spar',help='SPAR file',type=str)
        parser_philips.add_argument("-f", "--fileout", type=str,help="Output file base name (default = SDAT/SPAR name)")
        parser_philips.add_argument("-o", "--outdir", type=str,help="Output location (default = .)",default='.')
        parser_philips.add_argument('-j','--json',help='Create json sidecar.',action='store_true')
        parser_philips.set_defaults(func=self.philips)

        # Handle GE subcommand
        parser_ge = subparsers.add_parser('ge', help='Convert from GE p-file format.')
        parser_ge.add_argument('file',help='file to convert',type=str)
        parser_ge.add_argument("-f", "--fileout", type=str,help="Output file name (default = input name)")
        parser_ge.add_argument("-o", "--outdir", type=str,help="Output location (default = .)",default='.')
        parser_ge.add_argument('-j','--json',help='Create json sidecar.',action='store_true')
        parser_ge.set_defaults(func=self.ge)

        # Handle ismrmrd subcommand
        parser_ismrmrd = subparsers.add_parser('ismrmrd', help='Convert from ismrmrd format.')
        parser_ismrmrd.add_argument('file',help='file to convert',type=str)
        parser_ismrmrd.add_argument('a', type=str, help='placeholder')
        parser_ismrmrd.set_defaults(func=self.ismrmrd)

        # Handle text subcommand
        parser_txt = subparsers.add_parser('text', help='Convert from plain text format.')
        parser_txt.add_argument('file',help='file to convert',type=str)
        parser_txt.add_argument('-j','--json',help='Create json sidecar.',action='store_true')
        parser_txt.add_argument("-i", "--imagingfreq", type=float,help="Imaging (central) frequency in MHz",required=True)
        parser_txt.add_argument("-b", "--bandwidth", type=float,help="Reciever bandwidth (sweepwidth) in Hz.",required=True)
        parser_txt.add_argument("-a", "--affine", type=str,help="NIfTI affine file",required=False,metavar='<file>')
        parser_txt.add_argument("-f", "--fileout", type=str,help="Output file base name (default = input file name)")
        parser_txt.add_argument("-o", "--outdir", type=str,help="Output location (default = .)",default='.')        
        parser_txt.set_defaults(func=self.text)

        parser_jmrui = subparsers.add_parser('jmrui', help='Convert from jMRUI text format.')
        parser_jmrui.add_argument('file',help='file to convert',type=str)
        parser_jmrui.add_argument('-j','--json',help='Create json sidecar.',action='store_true')
        parser_jmrui.add_argument("-a", "--affine", type=str,help="NIfTI affine file",required=False,metavar='<file>')
        parser_jmrui.add_argument("-f", "--fileout", type=str,help="Output file base name (default = input file name)")
        parser_jmrui.add_argument("-o", "--outdir", type=str,help="Output location (default = .)",default='.')        
        parser_jmrui.set_defaults(func=self.jmrui)

        parser_raw = subparsers.add_parser('raw', help='Convert from LCModel RAW text format.')
        parser_raw.add_argument('file',help='file to convert',type=str)
        parser_raw.add_argument('-j','--json',help='Create json sidecar.',action='store_true')
        parser_raw.add_argument("-a", "--affine", type=str,help="NIfTI affine file",required=False,metavar='<file>')
        parser_raw.add_argument("-f", "--fileout", type=str,help="Output file base name (default = input file name)")
        parser_raw.add_argument("-o", "--outdir", type=str,help="Output location (default = .)",default='.')        
        parser_raw.set_defaults(func=self.raw)

        parser.add_argument('--verbose',action='store_true')        

        if len(sys.argv)==1:
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
            # Write nifti files            
            Path(self.outputDir).mkdir(parents=True, exist_ok=True)

            for n,i,o,d in zip(self.fileoutNames,self.imageOut,self.orientationInfoOut,self.dwellTimes):
                writeNii(n,self.outputDir,i,o,d)
            
            if self.metaData and args.json:
                for n,m in zip(self.fileoutNames,self.metaData):
                    writeJSON(n,self.outputDir,m)

        elif hasattr(args, 'view') and not args.view:
            print('No files to write.')

    def twix(self,args):     
        # Call mapVBVD to load the twix file.
        from mapVBVD import mapVBVD
        from spec2nii.twixfunctions import twix2DCMOrientation,examineTwix,extractTwixMetadata
        twixObj = mapVBVD(args.file,quiet=args.quiet)

        if  args.view:
            examineTwix(twixObj,op.basename(args.file),args.multiraid)
            return
        
        if isinstance(twixObj,list):                       
            twixObj = twixObj[args.multiraid-1]

        if not args.quiet:    
            print(f"Converting twix file {args.file}.")
            print(f'Looking for evalinfo flag {args.evalinfo}.')
        dataKey = args.evalinfo

        # Set squeeze data
        twixObj[dataKey].squeeze = True
        squeezedData = twixObj[dataKey]['']
        if not args.quiet: 
            print(f'Found data of size {squeezedData.shape}.')

        # Orientation calculations
        #1) Calculate dicom like imageOrientationPatient,imagePositionPatient,pixelSpacing and slicethickness
        imageOrientationPatient,imagePositionPatient,pixelSpacing,slicethickness = twix2DCMOrientation(twixObj['hdr'],verbose=args.verbose)
        # print(imageOrientationPatient)
        # print(imagePositionPatient)
        # print(pixelSpacing)
        # print(slicethickness)
        # 2) in style of dcm2niix
        # a) calculate Q44
        xyzMM = np.append(pixelSpacing,slicethickness)
        Q44 = nifti_dicom2mat(imageOrientationPatient,imagePositionPatient,xyzMM,verbose=args.verbose)
        # b) calculate nifti quaternion parameters
        Q44[:2,:] *= -1
        # qb,qc,qd,qx,qy,qz,dx,dy,dz,qfac = nifti_mat44_to_quatern(Q44)
        # 3) place in data class for nifti orientation parameters  
        currNiftiOrientation = NIFTIOrient(Q44)

        # Extract dwellTime
        dwellTime = twixObj['hdr']['MeasYaps'][('sRXSPEC','alDwellTime','0')]/1E9

        if args.json:
            # Extract metadata
            meta = extractTwixMetadata(twixObj['hdr'])

        # Identify what those indicies are
        # If cha is one: loop over 3rd and higher dims and make 2D images
        # If cha isn't present one: loop over 2nd and higher dims and make 1D images
        # Don't write here, just fill up class property lists for later writing
        if args.fileout:
            mainStr = args.fileout
        else:
            mainStr = op.splitext(op.basename(args.file))[0]

        dims = twixObj[dataKey].sqzDims()
        if dims[0] != 'Col':
            raise ValueError('Col is expected to be the first dimension in the Twix file, it is not.')# This is very unlikely to occur  but would cause complete failure.
        if 'Cha' in dims:
            if dims[1] != 'Cha':
                raise ValueError('If present Cha is expected to be the first dimension in the Twix file, it is not.') # This is very unlikely to occur but would cause complete failure.           
            if len(dims) ==2: # Only single file output, force singleton third dimensions
                squeezedData = squeezedData.reshape((squeezedData.shape+(1,)))
            # Loop over all the other dims from the third up
            for index in np.ndindex(squeezedData.shape[2:]):
                modIndex = (slice(None),slice(None),)+ index

                newshape = (1,1,1)+squeezedData[modIndex].shape
                data_out = squeezedData[modIndex].reshape(newshape)

                self.imageOut.append(data_out)
                self.orientationInfoOut.append(currNiftiOrientation)
                self.dwellTimes.append(dwellTime)
                if args.json:
                    self.metaData.append(meta)
                # Create strings
                for idx,ii in enumerate(index):
                    if len(dims) > 2:
                        indexStr = dims[2+idx]
                    else:
                        self.fileoutNames.append(f'{mainStr}')
                        continue

                    if idx==0:
                        self.fileoutNames.append(f'{mainStr}_{indexStr}{ii :03d}')
                    else:
                        self.fileoutNames[len(self.fileoutNames)-1] += f'_{indexStr}{ii :03d}'

        else: # Loop over all the other dims from the second up
            if len(dims) ==1: # Only single file output, force singleton second dimensions
                squeezedData = squeezedData.reshape((squeezedData.shape+(1,)))
            
            for index in np.ndindex(squeezedData.shape[1:]):
                modIndex = (slice(None),)+ index

                newshape = (1,1,1)+squeezedData[modIndex].shape
                data_out = squeezedData[modIndex].reshape(newshape)

                self.imageOut.append(data_out)
                self.orientationInfoOut.append(currNiftiOrientation)
                self.dwellTimes.append(dwellTime)
                if args.json:
                    self.metaData.append(meta)
                # Create strings
                for idx,ii in enumerate(index):
                    indexStr = dims[1+idx]
                    if idx==0:
                        self.fileoutNames.append(f'{mainStr}_{indexStr}{ii :03d}')
                    else:
                        self.fileoutNames[len(self.fileoutNames)-1] += f'_{indexStr}{ii :03d}'

    def dicom(self,args):
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


    def ismrmrd(self,args):
        print(f'Ismrmrd not yet handled!')
        print('exiting')
        return None

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