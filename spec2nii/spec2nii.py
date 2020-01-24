import argparse
import numpy as np
from dataclasses import dataclass
from spec2nii.writeNii import writeNii
from spec2nii.writeJSON import writeJSON
import os.path as op
from os import walk
from spec2nii.dcm2niiOrientation.orientationFuncs import nifti_dicom2mat,nifti_mat44_to_quatern

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
    Q44: np.array    

class spec2nii:
    def __init__(self):
        parser = argparse.ArgumentParser(description='Convert raw spectroscopy data to NIfTI format.')      

        subparsers = parser.add_subparsers(title='subcommands',
                                            description='File types supported')

        # Handle twix subcommand
        parser_twix = subparsers.add_parser('twix', help='Convert from Siemens .dat twix format.')
        parser_twix.add_argument('file',help='file to convert',type=str)
        parser_twix.add_argument('-j','--json',help='file to convert',action='store_true')
        group = parser_twix.add_mutually_exclusive_group(required=True)
        group.add_argument("-v", "--view", help="View contents of twix file, no files converted",action='store_true')
        group.add_argument('-e','--evalinfo', type=str, help='evalInfo flag to convert')
        parser_twix.add_argument("-f", "--fileout", type=str,help="Output file base name (default = input file name)")
        parser_twix.add_argument("-o", "--outdir", type=str,help="Output location (default = .)",default='.')
        parser_twix.add_argument("-m", "--multiraid", type=int,help="Select multiraid file to load (default = 2 i.e. 2nd file)",default=2)         
        parser_twix.set_defaults(func=self.twix)

        # Handle dicom subcommand
        parser_dicom = subparsers.add_parser('dicom', help='Convert from DICOM format.')
        parser_dicom.add_argument('file',help='file or directory to convert',type=str)
        parser_dicom.add_argument("-f", "--fileout", type=str,help="Output file base name (default = DCM SeriesDescription tag)")
        parser_dicom.add_argument("-o", "--outdir", type=str,help="Output location (default = .)",default='.')
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
        self.dwellTimes = []
        self.metaData = []

        self.outputDir = args.outdir

        args.func(args)

        if self.imageOut:
            # Write nifti files
            for n,i,o,d in zip(self.fileoutNames,self.imageOut,self.orientationInfoOut,self.dwellTimes):
                writeNii(n,self.outputDir,i,o,d)
            
            if self.metaData and args.json:
                for n,m in zip(self.fileoutNames,self.metaData):
                    writeJSON(n,self.outputDir,m)

        elif not args.view:
            print('No files to write.')

    def twix(self,args):     
        # Call mapVBVD to load the twix file.
        from mapVBVD import mapVBVD
        from spec2nii.twixfunctions import twix2DCMOrientation,examineTwix,extractTwixMetadata
        twixObj = mapVBVD.mapVBVD(self.fileIn)

        if  args.view:
            examineTwix(twixObj,op.basename(self.fileIn),args.multiraid)
            return
        
        if isinstance(twixObj,list):                       
            twixObj = twixObj[args.multiraid-1]
            
        print(f"Converting twix file {self.fileIn}.")
        print(f'Looking for evalinfo flag {args.evalinfo}.')
        dataKey = args.evalinfo

        # Set squeeze data
        twixObj[dataKey].squeeze = True
        squeezedData = twixObj[dataKey]['']
        print(f'Found data of size {squeezedData.shape}.')

        # Orientation calculations
        #1) Calculate dicom like imageOrientationPatient,imagePositionPatient,pixelSpacing and slicethickness
        imageOrientationPatient,imagePositionPatient,pixelSpacing,slicethickness = twix2DCMOrientation(twixObj['hdr'])
        # print(imageOrientationPatient)
        # print(imagePositionPatient)
        # print(pixelSpacing)
        # print(slicethickness)
        # 2) in style of dcm2niix
        # a) calculate Q44
        xyzMM = np.append(pixelSpacing,slicethickness)
        Q44 = nifti_dicom2mat(imageOrientationPatient,imagePositionPatient,xyzMM)
        # b) calculate nifti quaternion parameters
        Q44[:2,:] *= -1
        qb,qc,qd,qx,qy,qz,dx,dy,dz,qfac = nifti_mat44_to_quatern(Q44)
        # 3) place in data class for nifti orientation parameters  
        currNiftiOrientation = NIFTIOrient(qb,qc,qd,qx,qy,qz,dx,dy,dz,qfac,Q44)

        # Extract dwellTime
        dwellTime = twixObj['hdr']['MeasYaps'][('sRXSPEC','alDwellTime','0')]/1E9

        # Extract metadata
        meta = extractTwixMetadata(twixObj['hdr'])

        # Identify what those indicies are
        # If cha is one: loop over 3rd and higher dims and make 2D images
        # If cha isn't present one: loop over 2nd and higher dims and make 1D images
        # Don't write here, just fill up class property lists for later writing
        if args.fileout:
            mainStr = args.fileout
        else:
            mainStr = op.basename(self.fileIn)

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
                self.imageOut.append(squeezedData[modIndex])
                self.orientationInfoOut.append(currNiftiOrientation)
                self.dwellTimes.append(dwellTime)
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
                self.imageOut.append(squeezedData[modIndex])
                self.orientationInfoOut.append(currNiftiOrientation)
                self.dwellTimes.append(dwellTime)
                self.metaData.append(meta)
                # Create strings
                for idx,ii in enumerate(index):
                    indexStr = dims[1+idx]
                    if idx==0:
                        self.fileoutNames.append(f'{mainStr}_{indexStr}{ii :03d}')
                    else:
                        self.fileoutNames[len(self.fileoutNames)-1] += f'_{indexStr}{ii :03d}'

    def dicom(self,args):
        if op.isdir(self.fileIn):
            basePath = self.fileIn
            (_, _, filenames) = next(walk(basePath))
            print(f'Found {len(filenames)} files.')
        else:
            print('Single file conversion.')
            basePath = op.dirname(self.fileIn)
            filenames = [op.basename(self.fileIn),]
        
        # DICOM specific imports
        import nibabel.nicom.dicomwrappers
        for idx,fn in enumerate(filenames):
            print(f'Converting dicom file {self.fileIn}')
            
            img = nibabel.nicom.dicomwrappers.wrapper_from_file(op.join(basePath,fn))

            specData = np.frombuffer(img.dcm_data[('7fe1', '1010')].value, dtype=np.single)
            specDataCmplx = specData[0::2]+1j*specData[1::2]
            self.imageOut.append(specDataCmplx)

            if args.fileout:
                mainStr = args.fileout
            else: 
                mainStr = img.dcm_data.SeriesDescription

            self.fileoutNames.append(f'{mainStr}_{idx :03d}')

            #1) Extract dicom parameters
            imageOrientationPatient = np.array(img.csa_header['tags']['ImageOrientationPatient']['items']).reshape(2,3).transpose()
            imagePositionPatient = img.csa_header['tags']['VoiPosition']['items']
            xyzMM = np.array([img.csa_header['tags']['VoiPhaseFoV']['items'][0],
                         img.csa_header['tags']['VoiReadoutFoV']['items'][0],
                        img.csa_header['tags']['VoiThickness']['items'][0]])
            # 2) in style of dcm2niix
            # a) calculate Q44
            Q44 = nifti_dicom2mat(imageOrientationPatient,imagePositionPatient,xyzMM)
            # b) calculate nifti quaternion parameters
            Q44[:2,:] *= -1
            qb,qc,qd,qx,qy,qz,dx,dy,dz,qfac = nifti_mat44_to_quatern(Q44)
            # 3) place in data class for nifti orientation parameters 
            currNiftiOrientation = NIFTIOrient(qb,qc,qd,qx,qy,qz,dx,dy,dz,qfac,Q44)
            self.orientationInfoOut.append(currNiftiOrientation)
            self.dwellTimes.append(img.csa_header['tags']['RealDwellTime']['items'][0]*1E-9)
            self.metaData.append(extractDicomMetadata(img))


    def ismrmrd(self,args):
        print(f'Ismrmrd not yet handled!')
        print('exiting')
        return None

if __name__== "__main__":
    spec2nii()