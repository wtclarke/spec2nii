def extractDicomMetadata(dcmdata):
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
    metaDict.update({'SequenceName':dcmdata.csa_header['tags']['SequenceName']['items'][0]}})
    metaDict.update({'PulseSequenceDetails':''})
    metaDict.update({'EchoTime':dcmdata.csa_header['tags']['EchoTime']['items'][0]*1E-3})
    metaDict.update({'InversionTime':dcmdata.csa_header['tags']['InversionTime']['items'][0]})
    metaDict.update({'DwellTime':dcmdata.csa_header['tags']['RealDwellTime']['items'][0]*1E-9})
    metaDict.update({'FlipAngle':dcmdata.csa_header['tags']['FlipAngle']['items'][0]})
    metaDict.update({'InstitutionName':dcmdata.dcm_data.InstitutionName})
    metaDict.update({'InstitutionAddress':dcmdata.dcm_data.InstitutionAddress})
    metaDict.update({'InstitutionalDepartmentName':''})
    metaDict.update({'TimeStamp':dcmdata.dcm_data.AcquisitionTime})
    metaDict.update({'Nucleus':dcmdata.csa_header['tags']['ImagedNucleus']['items'][0]})
    metaDict.update({'ImagingFrequency':dcmdata.csa_header['tags']['ImagingFrequency']['items'][0]}})
    metaDict.update({'RepetitionTime':dcmdata.csa_header['tags']['RepetitionTime']['items'][0]/1E3})
    if dcmdata.csa_header['tags']['PatientOrientation']['items']:
        metaDict.update({'PatientPosition':dcmdata.csa_header['tags']['PatientOrientation']['items'][0]})
    metaDict.update({'ProtocolName':dcmdata.dcm_data.ProtocolName})
    metaDict.update({'ConversionSoftware':'spec2nii dicom'})
    return metaDict