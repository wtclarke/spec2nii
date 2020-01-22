import argparse
import numpy as np
import dataclasses
import writeNii
import os
# There are case specific imports below
@dataclass
class NIFTIOrient:
    qb: float
    qc: float
    qd: float
    qx: float
    qy: float
    qz: float
    dx: float
    dy: float
    dz: float
    qfac: float    

class spec2nii:
    def __init__(self):
        parser = argparse.ArgumentParser(description='Convert raw spectroscopy data to NIfTI format.')      

        subparsers = parser.add_subparsers(title='subcommands',
                                            description='File types supported')

        # Handle twix subcommand
        parser_twix = subparsers.add_parser('twix', help='Convert from Siemens .dat twix format.')
        parser_twix.add_argument('file',help='file to convert',type=str)
        parser_twix.add_argument('evalinfo', type=str, help='evalInfo flag to convert')
        parser_twix.add_argument("-f", "--fileout", type=str,help="Output file base name (default = input file name)")
        parser_twix.add_argument("-o", "--outdir", type=str,help="Output location (default = .)",default='.')
        parser_twix.set_defaults(func=self.twix)

        # Handle dicom subcommand
        parser_dicom = subparsers.add_parser('dicom', help='Convert from DICOM format.')
        parser_dicom.add_argument('file',help='file to convert',type=str)
        parser_dicom.add_argument('a', type=str, help='placeholder')
        parser_dicom.set_defaults(func=self.dicom)

        # Handle ismrmrd subcommand
        parser_ismrmrd = subparsers.add_parser('ismrmrd', help='Convert from ismrmrd format.')
        parser_ismrmrd.add_argument('file',help='file to convert',type=str)
        parser_ismrmrd.add_argument('a', type=str, help='placeholder')
        parser_ismrmrd.set_defaults(func=self.ismrmrd)

        args = parser.parse_args()
        
        self.fileIn = args.file
        self.fileoutNames = []
        self.imageOut = []
        self.orientationInfoOut = []

        self.outputDir = args.outdir

        args.func(args)

        # Write nifit files
        for n,i,o in zip(self.fileoutNames,self.imageOut,self.orientationInfoOut):
            writeNii(n,i,o)

    def twix(self,args):
        print(f"Converting twix file {self.fileIn}.")

        # Call mapVBVD to load the twix file.
        from mapVBVD import mapVBVD
        from twixfunctions import twix2DCMOrientation
        twixObj = mapVBVD.mapVBVD(self.fileIn)

        print(f'Looking for eval info flag {args.evalinfo}.')
        dataKey = args.evalinfo

        # Set squeeze data
        twixObj[dataKey].squeeze = True
        squeezedData = twixObj[dataKey]['']
        print(squeezedData.shape)

        # Orientation calculations
        #1) Calculate dicom like imageOrientationPatient,imagePositionPatient,pixelSpacing and slicethickness
        imageOrientationPatient,imagePositionPatient,pixelSpacing,slicethickness = twix2DCMOrientation(twixObj['hdr'])
        print(imageOrientationPatient)
        print(imagePositionPatient)
        print(pixelSpacing)
        print(slicethickness)
        # 2) in style of dcm2niix
        # a) calculate Q44
        xyzMM = np.append(pixelSpacing,slicethickness)
        Q44 = nifti_dicom2mat(imageOrientationPatient,imagePositionPatient,xyzMM)
        # b) calculate nifti quaternion parameters
        Q44[:2,:] *= -1
        qb,qc,qd,qx,qy,qz,dx,dy,dz,qfac = nifti_mat44_to_quatern(Q44)
        # 3) place in data class for nifti orientation parameters  
        currNiftiOrientation = NIFTIOrient(qb,qc,qd,qx,qy,qz,dx,dy,dz,qfac)

        # Identify what those indicies are
        # If cha is one: loop over 3rd and higher dims and make 2D images
        # If cha isn't present one: loop over 2nd and higher dims and make 1D images
        # Don't write here, just fill up class property lists for later writing
        if args.fileout:
            mainStr = args.fileout
        else:
            mainStr = os.path.basename(self.fileIn)
        print(mainStr)

        dims = twixObj['image'].sqzDims()
        if dims[0] != 'Col':
            raise ValueError('Col is expected to be the first dimension in the Twix file, it is not.')# This is very unlikely to occur
        if 'Cha' in dims:
            if dims[1] != 'Cha':
                raise ValueError('If present Cha is expected to be the first dimension in the Twix file, it is not.') # This is very unlikely to occur
            # Loop over all the other dims from the third up
            for index in np.ndindex(squeezedData.shape[2:]):
                modIndex = (slice(None),slice(None),)+ index
                self.imageOut.append(squeezedData[modIndex])
                self.orientationInfoOut.append(currNiftiOrientation)
                # Create strings
                for idx,ii in enumerate(index):
                    indexStr = dims[2+idx]
                    if idx==0:
                        self.fileoutNames.append(f'{mainStr}_{indexStr}{ii :03d}')
                    else:
                        self.fileoutNames[len(self.fileoutNames)-1] += f'_{indexStr}{ii :03d}'

        else: # Loop over all the other dims from the second up
            for index in np.ndindex(squeezedData.shape[1:]):
                modIndex = (slice(None),)+ index
                self.imageOut.append(squeezedData[modIndex])
                self.orientationInfoOut.append(currNiftiOrientation)
                # Create strings
                for idx,ii in enumerate(index):
                    indexStr = dims[1+idx]
                    if idx==0:
                        self.fileoutNames.append(f'{mainStr}_{indexStr}{ii :03d}')
                    else:
                        self.fileoutNames[len(self.fileoutNames)-1] += f'_{indexStr}{ii :03d}'

    def dicom(self,args):
        print(f'Converting dicom file {self.fileIn}')

    def ismrmrd(self,args):
        print(f'Converting ismrmrd file {self.fileIn}')

if __name__== "__main__":
    spec2nii()