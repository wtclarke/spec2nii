"""spec2nii module containing functions specific to interpreting Philips DICOM
Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford
"""
import numpy as np
import nibabel.nicom.dicomwrappers
from spec2nii.dcm2niiOrientation.orientationFuncs import dcm_to_nifti_orientation
from spec2nii.nifti_orientation import NIFTIOrient
from spec2nii import nifti_mrs
from datetime import datetime


class CSINotHandledError(Exception):
    pass


def svs_or_CSI(img):
    """Identify from the headers whether data is CSI, SVS or FID."""
    if img.image_shape is None:
        if img.dcm_data.VolumeLocalizationTechnique == 'NONE':
            return 'FID'
        else:
            return 'SVS'
    elif np.product(img.image_shape) > 1.0:
        return 'CSI'
    else:
        return 'SVS'


def multi_file_dicom(files_in, fname_out, tag, verbose):
    """Parse a list of UIH DICOM files"""

    # Convert each file (combine after)
    data_list = []
    ref_list = []
    orientation_list = []
    dwelltime_list = []
    meta_list = []
    mainStr = ''
    for idx, fn in enumerate(files_in):
        if verbose:
            print(f'Converting dicom file {fn}')

        img = nibabel.nicom.dicomwrappers.wrapper_from_file(fn)

        mrs_type = svs_or_CSI(img)

        if mrs_type == 'SVS':
            specDataCmplx, spec_ref, orientation, dwelltime, meta_obj = _process_philips_svs(img, verbose)

            newshape = (1, 1, 1) + specDataCmplx.shape
            spec_data = specDataCmplx.reshape(newshape)
            spec_ref = spec_ref.reshape(newshape)
        elif mrs_type == 'FID':
            specDataCmplx, spec_ref, orientation, dwelltime, meta_obj = _process_philips_fid(img, verbose)

            newshape = (1, 1, 1) + specDataCmplx.shape
            spec_data = specDataCmplx.reshape(newshape)
            spec_ref = spec_ref.reshape(newshape)
        else:
            raise CSINotHandledError('CSI data is currently not handled for the Philips DICOM format.'
                                     'Please contact the developers if you have examples of this type of data.')

        data_list.append(spec_data)
        ref_list.append(spec_ref)
        orientation_list.append(orientation)
        dwelltime_list.append(dwelltime)
        meta_list.append(meta_obj)

        if idx == 0:
            if fname_out:
                mainStr = fname_out
            else:
                mainStr = img.dcm_data.SeriesDescription

    # If data shape, orientation and dwelltime match combine
    # into one NIFTI MRS object.
    # Otherwise return a list of files/names
    def all_equal(lst):
        return lst[:-1] == lst[1:]

    combine = all_equal([d.shape for d in data_list])\
        and all_equal([o.Q44.tolist() for o in orientation_list])\
        and all_equal(dwelltime_list)

    nifti_mrs_out, fnames_out = [], []
    if combine:
        # Combine files into single MRS NIfTI
        # Single file name
        fnames_out.append(mainStr)
        fnames_out.append(mainStr + '_ref')

        dt_used = dwelltime_list[0]
        or_used = orientation_list[0]

        # Add original files to nifti meta information.
        meta_used = meta_list[0]
        meta_used.set_standard_def('OriginalFile', [str(ff.name) for ff in files_in])

        # Combine data into 5th dimension if needed
        if len(data_list) > 1:
            combined_data = np.stack(data_list, axis=-1)
            combined_ref = np.stack(ref_list, axis=-1)
        else:
            combined_data = data_list[0]
            combined_ref = ref_list[0]

        # Add dimension information (if not None for default)
        if tag:
            meta_used.set_dim_info(0, tag)

        # Create NIFTI MRS object.
        nifti_mrs_out.append(nifti_mrs.NIfTI_MRS(combined_data, or_used.Q44, dt_used, meta_used))
        nifti_mrs_out.append(nifti_mrs.NIfTI_MRS(combined_ref, or_used.Q44, dt_used, meta_used))
    else:
        for idx, (dd, rr, oo, dt, mm, ff) in enumerate(zip(data_list,
                                                       ref_list,
                                                       orientation_list,
                                                       dwelltime_list,
                                                       meta_list,
                                                       files_in)):
            # Add original files to nifti meta information.
            mm.set_standard_def('OriginalFile', [str(ff.name), ])
            fnames_out.append(f'{mainStr}_{idx:03}')
            nifti_mrs_out.append(nifti_mrs.NIfTI_MRS(dd, oo.Q44, dt, mm))

            fnames_out.append(f'{mainStr}_ref_{idx:03}')
            nifti_mrs_out.append(nifti_mrs.NIfTI_MRS(rr, oo.Q44, dt, mm))

    return nifti_mrs_out, fnames_out


def _process_philips_svs(img, verbose):
    """Process Philips DICOM SVS data"""

    specData = np.frombuffer(img.dcm_data[('5600', '0020')].value, dtype=np.single)
    specDataCmplx = specData[0::2] + 1j * specData[1::2]

    # In the one piece of data I have been provided the data is twice as long as indicated (1 avg)
    # and the second half is a water reference.
    spec_points = img.dcm_data.SpectroscopyAcquisitionDataColumns
    spec_data_main = specDataCmplx[:spec_points]
    spec_data_ref = specDataCmplx[spec_points:]

    # 1) Extract dicom parameters
    currNiftiOrientation = _enhanced_dcm_svs_to_orientation(img, verbose)

    dwelltime = 1.0 / img.dcm_data.SpectralWidth
    meta = _extractDicomMetadata(img)

    return spec_data_main, spec_data_ref, currNiftiOrientation, dwelltime, meta


def _process_philips_fid(img, verbose):
    """Process Philips DICOM FID data"""

    specData = np.frombuffer(img.dcm_data[('5600', '0020')].value, dtype=np.single)
    specDataCmplx = specData[0::2] + 1j * specData[1::2]

    # In the one piece of data I have been provided the data is twice as long as indicated (1 avg)
    # and the second half is a water reference.
    spec_points = img.dcm_data.SpectroscopyAcquisitionDataColumns
    spec_data_main = specDataCmplx[:spec_points]
    spec_data_ref = specDataCmplx[spec_points:]

    # 1) Extract dicom parameters
    defaultaffine = np.diag(np.array([10000, 10000, 10000, 1]))
    currNiftiOrientation = NIFTIOrient(defaultaffine)

    dwelltime = 1.0 / img.dcm_data.SpectralWidth
    meta = _extractDicomMetadata(img)

    return spec_data_main, spec_data_ref, currNiftiOrientation, dwelltime, meta


def _enhanced_dcm_svs_to_orientation(img, verbose):
    '''Convert the Volume Localization Sequence (0018,9126) enhanced DICOM tag
    to a 4x4 affine format.

    Assumes the center of all slabs are aligned currently.'''

    # affine3 = []
    # voxsize = []
    # for vl in img.dcm_data.VolumeLocalizationSequence:
    #     affine3.append(vl.SlabOrientation)
    #     voxsize.append(vl.SlabThickness)
    # affine3 = np.asarray(affine3)
    # voxsize = np.diag(np.asarray(voxsize))

    # locSeq = img.dcm_data.VolumeLocalizationSequence
    # if not np.allclose(locSeq[0].MidSlabPosition, locSeq[1].MidSlabPosition) \
    #         or not np.allclose(locSeq[0].MidSlabPosition, locSeq[2].MidSlabPosition):
    #     raise ValueError('Mid Slab Position must be the same for all three entries')

    # pos3 = np.asarray(img.dcm_data.VolumeLocalizationSequence[0].MidSlabPosition)

    # affine4 = np.eye(4)
    # affine4[:3, :3] = affine3 @ voxsize
    # affine4[:3, 3] = pos3

    # Above didn't work. Drop back to original method.

    # Extract dicom parameters
    imageOrientationPatient = img.dcm_data.PerFrameFunctionalGroupsSequence[0]\
                                 .PlaneOrientationSequence[0].ImageOrientationPatient
    imageOrientationPatient = np.asarray(imageOrientationPatient).reshape(2, 3)
    imagePositionPatient = img.dcm_data.PerFrameFunctionalGroupsSequence[0]\
                              .PlanePositionSequence[0].ImagePositionPatient
    imagePositionPatient = np.asarray(imagePositionPatient)
    xyzMM = [float(img.dcm_data.PerFrameFunctionalGroupsSequence[0].PixelMeasuresSequence[0].PixelSpacing[0]),
             float(img.dcm_data.PerFrameFunctionalGroupsSequence[0].PixelMeasuresSequence[0].PixelSpacing[1]),
             float(img.dcm_data.PerFrameFunctionalGroupsSequence[0].PixelMeasuresSequence[0].SliceThickness)]
    xyzMM = np.asarray(xyzMM)

    currNiftiOrientation = dcm_to_nifti_orientation(imageOrientationPatient,
                                                    imagePositionPatient,
                                                    xyzMM,
                                                    (1, 1, 1),
                                                    half_shift=False,
                                                    verbose=verbose)

    return currNiftiOrientation


def _extractDicomMetadata(dcmdata):
    """ Extract information from the nibabel DICOM object to insert into the json header ext.

    There seems to be a large 'uprotocol' in tag 0065,1007 but I don't know how to interpret it.

    Args:
        dcmdata: nibabel.nicom image object
    Returns:
        obj (hdr_ext): NIfTI MRS hdr ext object.
    """

    # Extract required metadata and create hdr_ext object
    obj = nifti_mrs.hdr_ext(dcmdata.dcm_data.TransmitterFrequency,
                            dcmdata.dcm_data.ResonantNucleus)

    # Some scanner information
    obj.set_standard_def('Manufacturer', dcmdata.dcm_data.Manufacturer)
    obj.set_standard_def('ManufacturersModelName', dcmdata.dcm_data.ManufacturerModelName)
    if 'DeviceSerialNumber' in dcmdata.dcm_data:
        obj.set_standard_def('DeviceSerialNumber', str(dcmdata.dcm_data.DeviceSerialNumber))
    obj.set_standard_def('SoftwareVersions', str(dcmdata.dcm_data.SoftwareVersions))

    obj.set_standard_def('InstitutionName', dcmdata.dcm_data.InstitutionName)
    if 'InstitutionAddress' in dcmdata.dcm_data:
        obj.set_standard_def('InstitutionAddress', dcmdata.dcm_data.InstitutionAddress)

    # Rx Coil name in header.
    # ToDo
    # img.dcm_data.SharedFunctionalGroupsSequence[0].MRReceiveCoilSequence[0].ReceiveCoilName

    # Some sequence information
    obj.set_standard_def('SequenceName', dcmdata.dcm_data.PulseSequenceName)
    # obj.set_standard_def('ProtocolName', dcmdata.dcm_data.ProtocolName)

    # Some subject information
    obj.set_standard_def('PatientPosition', dcmdata.dcm_data.PatientPosition)

    def set_if_present(tag, val):
        if val:
            obj.set_standard_def(tag, val)
    set_if_present('PatientDoB', dcmdata.dcm_data.PatientBirthDate)
    set_if_present('PatientSex', dcmdata.dcm_data.PatientSex)

    if dcmdata.dcm_data.PatientName:
        obj.set_standard_def('PatientName', dcmdata.dcm_data.PatientName.family_name)

    if 'PatientWeight' in dcmdata.dcm_data and dcmdata.dcm_data.PatientWeight:
        obj.set_standard_def('PatientWeight', float(dcmdata.dcm_data.PatientWeight))

    # Timing and sequence parameters
    echo_time = float(dcmdata.dcm_data.PerFrameFunctionalGroupsSequence[0]
                      .MREchoSequence[0].EffectiveEchoTime) * 1E-3
    rep_tim = float(dcmdata.dcm_data.PerFrameFunctionalGroupsSequence[0]
                    .MRTimingAndRelatedParametersSequence[0].RepetitionTime) / 1E3
    fa = float(dcmdata.dcm_data.PerFrameFunctionalGroupsSequence[0]
               .MRTimingAndRelatedParametersSequence[0].RepetitionTime)
    obj.set_standard_def('EchoTime', echo_time)
    obj.set_standard_def('ExcitationFlipAngle', fa)
    obj.set_standard_def('RepetitionTime', rep_tim)
    # TO DO  - nibabel might need updating.
    # obj.set_standard_def('TxOffset', )

    # Conversion information
    obj.set_standard_def('ConversionMethod', 'spec2nii')
    conversion_time = datetime.now().isoformat(sep='T', timespec='milliseconds')
    obj.set_standard_def('ConversionTime', conversion_time)

    return obj
