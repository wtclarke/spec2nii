""" spec2nii module containing functions specific to interpreting Siemens DICOM
Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford 
"""
import numpy as np
from spec2nii.dcm2niiOrientation.orientationFuncs import nifti_dicom2mat
from spec2nii.nifti_orientation import NIFTIOrient

def svs_or_CSI(img):
    rows = img.csa_header['tags']['Rows']['items'][0]
    cols = img.csa_header['tags']['Columns']['items'][0]
    slices = img.csa_header['tags']['NumberOfFrames']['items'][0]

    if np.prod([rows,cols,slices])>1.0:
        return 'CSI'
    else:
        return 'SVS'

def process_siemens_svs(img,args):
    specData = np.frombuffer(img.dcm_data[('7fe1', '1010')].value, dtype=np.single)
    specDataCmplx = specData[0::2]+1j*specData[1::2]    

    #1) Extract dicom parameters
    imageOrientationPatient = np.array(img.csa_header['tags']['ImageOrientationPatient']['items']).reshape(2,3)
    imagePositionPatient = img.csa_header['tags']['VoiPosition']['items']
    xyzMM = np.array([img.csa_header['tags']['VoiPhaseFoV']['items'][0],
                    img.csa_header['tags']['VoiReadoutFoV']['items'][0],
                img.csa_header['tags']['VoiThickness']['items'][0]])
    # 2) in style of dcm2niix
    # a) calculate Q44
    Q44 = nifti_dicom2mat(imageOrientationPatient,imagePositionPatient,xyzMM)
    # b) calculate nifti quaternion parameters
    Q44[:2,:] *= -1
    # qb,qc,qd,qx,qy,qz,dx,dy,dz,qfac = nifti_mat44_to_quatern(Q44)
    # 3) place in data class for nifti orientation parameters 
    currNiftiOrientation = NIFTIOrient(Q44)
    dwelltime = img.csa_header['tags']['RealDwellTime']['items'][0]*1E-9
    meta = extractDicomMetadata(img)

    return specDataCmplx,currNiftiOrientation,dwelltime,meta

def process_siemens_csi(img,args):
    specData = np.frombuffer(img.dcm_data[('7fe1', '1010')].value, dtype=np.single)
    specDataCmplx = specData[0::2]-1j*specData[1::2]

    rows = img.csa_header['tags']['Rows']['items'][0]
    cols = img.csa_header['tags']['Columns']['items'][0]
    slices = img.csa_header['tags']['NumberOfFrames']['items'][0]
    spectral_points = img.csa_header['tags']['DataPointColumns']['items'][0]    

    specDataCmplx = specDataCmplx.reshape((slices,rows,cols,spectral_points))
    specDataCmplx = np.moveaxis(specDataCmplx,(0,1,2),(2,1,0))

    #1) Extract dicom parameters
    imageOrientationPatient = np.array(img.csa_header['tags']['ImageOrientationPatient']['items']).reshape(2,3)
    imagePositionPatient = np.array(img.csa_header['tags']['ImagePositionPatient']['items'])
    xyzMM = np.array([img.csa_header['tags']['PixelSpacing']['items'][0],
                  img.csa_header['tags']['PixelSpacing']['items'][1],
                  img.csa_header['tags']['SliceThickness']['items'][0]])
    # 2) in style of dcm2niix
    # a) calculate Q44
    Q44 = nifti_dicom2mat(imageOrientationPatient,imagePositionPatient,xyzMM)
    # b) calculate nifti quaternion parameters
    Q44[:2,:] *= -1
    # qb,qc,qd,qx,qy,qz,dx,dy,dz,qfac = nifti_mat44_to_quatern(Q44)
    # 3) place in data class for nifti orientation parameters 
    currNiftiOrientation = NIFTIOrient(Q44)
    dwelltime = img.csa_header['tags']['RealDwellTime']['items'][0]*1E-9
    meta = extractDicomMetadata(img)
    
    return specDataCmplx,currNiftiOrientation,dwelltime,meta


def extractDicomMetadata(dcmdata):
    """ Extract information from the nibabel DICOM objhect to insert into the json sidecar.
    
    Args:
        dcmdata: nibabel.nicom image object
    Returns:
        metaDict (dict): Json sidecard output
    """

    metaDict = {}
    metaDict.update({'Modality':'MR'})
    metaDict.update({'Manufacturer':dcmdata.dcm_data.Manufacturer})
    metaDict.update({'ManufacturersModelName':dcmdata.dcm_data.ManufacturerModelName})
    metaDict.update({'DeviceSerialNumber':dcmdata.dcm_data.DeviceSerialNumber})
    metaDict.update({'StationName':dcmdata.dcm_data.StationName}) # Not an exact match but close
    metaDict.update({'SoftwareVersions':dcmdata.dcm_data.SoftwareVersions})
    metaDict.update({'MagneticFieldStrength':dcmdata.csa_header['tags']['MagneticFieldStrength']['items'][0]})
    if len(dcmdata.csa_header['tags']['ReceivingCoil']['items'])>0:
        metaDict.update({'ReceiveCoilName':dcmdata.csa_header['tags']['ReceivingCoil']['items'][0]})
    else:
         metaDict.update({'ReceiveCoilName':dcmdata.csa_header['tags']['ImaCoilString']['items'][0]})
    
    metaDict.update({'ScanningSequence':dcmdata.csa_header['tags']['ScanningSequence']['items'][0]})
    metaDict.update({'SequenceVariant':''})
    metaDict.update({'ScanOptions':''})
    metaDict.update({'SequenceName':dcmdata.csa_header['tags']['SequenceName']['items'][0]})
    metaDict.update({'PulseSequenceDetails':''})
    metaDict.update({'EchoTime':dcmdata.csa_header['tags']['EchoTime']['items'][0]*1E-3})
    try: 
        metaDict.update({'InversionTime':dcmdata.csa_header['tags']['InversionTime']['items'][0]})
    except:
        metaDict.update({'InversionTime':''})
    metaDict.update({'Dwelltime':dcmdata.csa_header['tags']['RealDwellTime']['items'][0]*1E-9})
    metaDict.update({'FlipAngle':dcmdata.csa_header['tags']['FlipAngle']['items'][0]})
    metaDict.update({'InstitutionName':dcmdata.dcm_data.InstitutionName})
    metaDict.update({'InstitutionAddress':dcmdata.dcm_data.InstitutionAddress})
    metaDict.update({'InstitutionalDepartmentName':''})
    metaDict.update({'TimeStamp':dcmdata.dcm_data.AcquisitionTime})
    metaDict.update({'Nucleus':dcmdata.csa_header['tags']['ImagedNucleus']['items'][0]})
    metaDict.update({'ImagingFrequency':dcmdata.csa_header['tags']['ImagingFrequency']['items'][0]})
    metaDict.update({'RepetitionTime':dcmdata.csa_header['tags']['RepetitionTime']['items'][0]/1E3})
    if dcmdata.csa_header['tags']['PatientOrientation']['items']:
        metaDict.update({'PatientPosition':dcmdata.csa_header['tags']['PatientOrientation']['items'][0]})
    metaDict.update({'ProtocolName':dcmdata.dcm_data.ProtocolName})
    metaDict.update({'ConversionSoftware':'spec2nii dicom'})
    return metaDict