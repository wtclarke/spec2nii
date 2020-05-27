""" spec2nii module containing functions specific to Siemens TWIX format
Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford 
"""

import spec2nii.GSL.gslfunctions as GSL
import numpy as np

def twix2DCMOrientation(mapVBVDHdr,verbose=False):
    """ Convert twix orientation information to DICOM equivalent.
    
    Convert orientation to DICOM imageOrientationPatient, imagePositionPatient, 
    pixelSpacing and sliceThickness field values.

    Args:
        mapVBVDHdr (dict): Header info interpreted by pymapVBVD
        verbose (bool,optionl)
    Returns:
        imageOrientationPatient
        imagePositionPatient
        pixelSpacing
        sliceThickness

    """
    #Orientation information
    if ('sSpecPara','sVoI','sNormal','dSag') in mapVBVDHdr['MeasYaps']:
        NormaldSag = mapVBVDHdr['MeasYaps'][('sSpecPara','sVoI','sNormal','dSag')]
    else: 
        NormaldSag = 0.0
    
    if ('sSpecPara','sVoI','sNormal','dCor') in mapVBVDHdr['MeasYaps']:
        NormaldCor = mapVBVDHdr['MeasYaps'][('sSpecPara','sVoI','sNormal','dCor')]
    else :
        NormaldCor = 0.0

    if ('sSpecPara','sVoI','sNormal','dTra') in mapVBVDHdr['MeasYaps']:
        NormaldTra = mapVBVDHdr['MeasYaps'][('sSpecPara','sVoI','sNormal','dTra')]
    else:
        NormaldTra = 0.0
    
    if ('sSpecPara','sVoI','dInPlaneRot') in mapVBVDHdr['MeasYaps']: 
        inplaneRotation = mapVBVDHdr['MeasYaps'][('sSpecPara','sVoI','dInPlaneRot')]
    else:
        inplaneRotation = 0.0    

    TwixSliceNormal =np.array([NormaldSag,NormaldCor,NormaldTra],dtype = float)

    RoFoV = mapVBVDHdr['MeasYaps'][('sSpecPara','sVoI','dReadoutFOV')]
    PeFoV = mapVBVDHdr['MeasYaps'][('sSpecPara','sVoI','dPhaseFOV')]              

    dColVec_vector, dRowVec_vector = GSL.fGSLCalcPRS(TwixSliceNormal,inplaneRotation,verbose)

    imageOrientationPatient = np.stack((dRowVec_vector,dColVec_vector),axis = 0)
    sliceNormal = TwixSliceNormal

    columns = 1
    rows = 1
    slices= 1

    pixelSpacing = np.array([PeFoV, RoFoV]) #[RoFoV PeFoV];
    sliceThickness = mapVBVDHdr['MeasYaps'][('sSpecPara','sVoI','dThickness')]

    # Position info
    if ('sSpecPara','sVoI','sPosition','dSag') in mapVBVDHdr['MeasYaps']:
        PosdSag = mapVBVDHdr['MeasYaps'][('sSpecPara','sVoI','sPosition','dSag')]
    else: 
        PosdSag = 0.0
    
    if ('sSpecPara','sVoI','sPosition','dCor') in mapVBVDHdr['MeasYaps']:
        PosdCor = mapVBVDHdr['MeasYaps'][('sSpecPara','sVoI','sPosition','dCor')]
    else:
        PosdCor = 0.0

    if ('sSpecPara','sVoI','sPosition','dTra') in mapVBVDHdr['MeasYaps']:
        PosdTra = mapVBVDHdr['MeasYaps'][('sSpecPara','sVoI','sPosition','dTra')]
    else:
        PosdTra = 0.0

    basePosition =np.array([PosdSag, PosdCor, PosdTra],dtype = float)
    imagePositionPatient = basePosition
    if verbose:
        print(f'imagePositionPatient is {imagePositionPatient.ravel()}')
        print(f'imageOrientationPatient is \n{imageOrientationPatient}')
        print(f'{imageOrientationPatient.ravel()}')
        print(f'pixelSpacing is {pixelSpacing}')
        
    return imageOrientationPatient,imagePositionPatient,pixelSpacing,sliceThickness

def examineTwix(twixObj,fileName,mraid):
    """ Print formated twix contents"""

    print(f'Contents of file: {fileName}')

    if isinstance(twixObj,list):
        print(f'Multiraid file, {len(twixObj)} files found.')
        print(f'Selecting file {mraid}. Use -m option to change.')                       
        twixObj = twixObj[mraid-1]

    evalInfoFlags = twixObj.keys()
    evalInfoFlags = [i for i in evalInfoFlags if i != 'hdr']

    print(f'The file contains these evalinfo flags with dimensions and sizes as follows:')
    for ev in evalInfoFlags:
        twixObj[ev].squeeze = True
        tmpSqzSize = twixObj[ev].sqzSize()
        tmpSqzDims = ', '.join(twixObj[ev].sqzDims())
        print(f'{ev: <15}:\t{tmpSqzDims: <20}\t{tmpSqzSize}')

    #import pdb; pdb.set_trace()

def extractTwixMetadata(mapVBVDHdr):
    """ Extract information from the pymapVBVD header to insert into the json sidecar.
    
    Args:
        dcmdata (dict): Twix headers
    Returns:
        metaDict (dict): Json sidecard output
    """

    metaDict = {}
    metaDict.update({'Modality':'MR'})
    metaDict.update({'Manufacturer':mapVBVDHdr['Dicom'][('Manufacturer')]})
    metaDict.update({'ManufacturersModelName':mapVBVDHdr['Dicom'][('ManufacturersModelName')]})
    metaDict.update({'DeviceSerialNumber':_try_int(mapVBVDHdr['Dicom'][('DeviceSerialNumber')])})
    metaDict.update({'StationName':_try_int(mapVBVDHdr['Dicom'][('DeviceSerialNumber')])}) # Not an exact match but close
    metaDict.update({'SoftwareVersions':mapVBVDHdr['Dicom'][('SoftwareVersions')]})
    metaDict.update({'MagneticFieldStrength':mapVBVDHdr['Dicom']['flMagneticFieldStrength']})
    if ('sCoilSelectMeas', 'aRxCoilSelectData', '0', 'asList', '0', 'sCoilElementID', 'tCoilID') in mapVBVDHdr['MeasYaps']:
        metaDict.update({'ReceiveCoilName':mapVBVDHdr['MeasYaps'][('sCoilSelectMeas', 'aRxCoilSelectData', '0', 'asList', '0', 'sCoilElementID', 'tCoilID')]})
    elif ('asCoilSelectMeas', '0', 'asList', '0', 'sCoilElementID', 'tCoilID') in  mapVBVDHdr['MeasYaps']:
         metaDict.update({'ReceiveCoilName':mapVBVDHdr['MeasYaps'][('asCoilSelectMeas', '0', 'asList', '0', 'sCoilElementID', 'tCoilID')]})
    else:
        metaDict.update({'ReceiveCoilName':''})
    metaDict.update({'ScanningSequence':mapVBVDHdr['Dicom'][('tScanningSequence')]})
    metaDict.update({'SequenceVariant':mapVBVDHdr['Dicom'][('tSequenceVariant')]})
    metaDict.update({'ScanOptions':mapVBVDHdr['Dicom'][('tScanOptions')]})
    metaDict.update({'SequenceName':''}) # Can't find in twix headers
    metaDict.update({'PulseSequenceDetails':mapVBVDHdr['Config'][('SequenceFileName')]})
    metaDict.update({'EchoTime':mapVBVDHdr['Phoenix'][('alTE','0')]*1E-6})
    if ('InversionTime') in mapVBVDHdr['Meas']:
        metaDict.update({'InversionTime':mapVBVDHdr['Meas'][('TI_Time')]})
    else:
        metaDict.update({'InversionTime':''})    
    metaDict.update({'Dwelltime':mapVBVDHdr['Meas'][('DwellTimeSig')]*1e-9})
    metaDict.update({'FlipAngle':mapVBVDHdr['Meas'][('FlipAngle')]})
    metaDict.update({'InstitutionName':mapVBVDHdr['Dicom'][('InstitutionName')]})
    metaDict.update({'InstitutionAddress':mapVBVDHdr['Dicom'][('InstitutionAddress')]})
    metaDict.update({'InstitutionalDepartmentName':''})#Not in twix headers
    metaDict.update({'TimeStamp':''})#Not in twix headers
    metaDict.update({'Nucleus':mapVBVDHdr['Meas'][('ResonantNucleus')]})
    metaDict.update({'ImagingFrequency':mapVBVDHdr['Meas'][('Frequency')]/1E6})
    if ('TR_Time') in mapVBVDHdr['Meas']:
        metaDict.update({'RepetitionTime':mapVBVDHdr['Meas'][('TR_Time')]/1E6}) 
    else:
        metaDict.update({'RepetitionTime':mapVBVDHdr['Meas'][('TR')]/1E6})     
    metaDict.update({'PatientPosition':mapVBVDHdr['Meas'][('PatientPosition')]})
    metaDict.update({'ProtocolName':mapVBVDHdr['Dicom'][('tProtocolName')]})
    metaDict.update({'ConversionSoftware':'spec2nii'})
    return metaDict

def _try_int(value):

    try:
        return int(value)
    except:
        return value