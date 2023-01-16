"""spec2nii module containing functions specific to interpreting Philips DICOM
Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford
"""
from datetime import datetime

import numpy as np
import nibabel.nicom.dicomwrappers

from nifti_mrs.create_nmrs import gen_nifti_mrs_hdr_ext
from nifti_mrs.hdr_ext import Hdr_Ext

from spec2nii.dcm2niiOrientation.orientationFuncs import dcm_to_nifti_orientation
from spec2nii.nifti_orientation import NIFTIOrient
from spec2nii import __version__ as spec2nii_ver


class CSINotHandledError(Exception):
    pass


class IncorrectDataDimensionsError(Exception):
    pass


def _is_new_format(img):
    """Is the file old or new format DICOM"""
    media_storage_class = img.dcm_data.file_meta[('0002', '0002')].value
    if media_storage_class == '1.2.840.10008.5.1.4.1.1.4.2':
        return True
    elif media_storage_class == '1.3.46.670589.11.0.0.12.1':
        return False
    else:
        raise ValueError(f'Media Storage SOP Class UID {media_storage_class} not recognised as spectroscopy.')


def svs_or_CSI(img):
    """Identify from the headers whether data is CSI, SVS or FID."""
    if _is_new_format(img):
        if img.image_shape is None \
                and img.dcm_data.VolumeLocalizationTechnique == 'NONE':
            return 'FID'
        elif np.product(img.image_shape) > 1.0:
            return 'CSI'
        else:
            return 'SVS'
    else:
        private_mrs_tags = img.dcm_data[('2005', '140f')][0]
        columns = private_mrs_tags[('0018', '9234')].value  # 'Spectroscopy Acquisition Phase Columns'
        rows = private_mrs_tags[('0018', '9095')].value  # 'Spectroscopy Acquisition Phase Rows'
        slices = private_mrs_tags[('0018', '9159')].value  # 'Spectroscopy Acquisition Phase Out-of-plane'

        if img.image_shape is None \
                and img.dcm_data.VolumeLocalizationTechnique == 'NONE':
            return 'FID'
        elif np.product([columns, rows, slices]) > 1.0:
            return 'CSI'
        else:
            return 'SVS'


def multi_file_dicom(files_in, fname_out, tag, verbose):
    """Parse a list of Philips DICOM files"""

    # Convert each file (combine after)
    data_list = []
    ref_list = []
    orientation_list = []
    dwelltime_list = []
    meta_list = []
    meta_ref_list = []
    mainStr = ''
    for idx, fn in enumerate(files_in):
        if verbose:
            print(f'Converting dicom file {fn}')

        img = nibabel.nicom.dicomwrappers.wrapper_from_file(fn)

        mrs_type = svs_or_CSI(img)

        if mrs_type == 'SVS':
            specDataCmplx, spec_ref, orientation, dwelltime, meta_obj, meta_ref_obj = _process_philips_svs(img, verbose)

            newshape = (1, 1, 1) + specDataCmplx.shape
            spec_data = specDataCmplx.reshape(newshape)
            if spec_ref is not None:
                spec_ref = spec_ref.reshape(newshape)

            # Data appears to require conjugation to meet standard's conventions.
            spec_data = spec_data.conj()
            if spec_ref is not None:
                spec_ref = spec_ref.conj()

        elif mrs_type == 'FID':
            specDataCmplx, spec_ref, orientation, dwelltime, meta_obj, meta_ref_obj = _process_philips_fid(img, verbose)

            newshape = (1, 1, 1) + specDataCmplx.shape
            spec_data = specDataCmplx.reshape(newshape)
            spec_ref = spec_ref.reshape(newshape)

            # Data appears to require conjugation to meet standard's conventions.
            spec_data = spec_data.conj()
            spec_ref = spec_ref.conj()
        else:
            raise CSINotHandledError('CSI data is currently not handled for the Philips DICOM format.'
                                     'Please contact the developers if you have examples of this type of data.')

        data_list.append(spec_data)
        ref_list.append(spec_ref)
        orientation_list.append(orientation)
        dwelltime_list.append(dwelltime)
        meta_list.append(meta_obj)
        meta_ref_list.append(meta_ref_obj)

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
        if ref_list[0] is not None:
            fnames_out.append(mainStr + '_ref')

        dt_used = dwelltime_list[0]
        or_used = orientation_list[0]

        # Add original files to nifti meta information.
        meta_used = meta_list[0]
        meta_used.set_standard_def('OriginalFile', [str(ff.name) for ff in files_in])
        if meta_ref_list[0] is not None:
            meta_ref_used = meta_ref_list[0]
            meta_ref_used.set_standard_def('OriginalFile', [str(ff.name) for ff in files_in])

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
        nifti_mrs_out.append(gen_nifti_mrs_hdr_ext(combined_data, dt_used, meta_used, or_used.Q44, no_conj=True))
        if ref_list[0] is not None:
            nifti_mrs_out.append(gen_nifti_mrs_hdr_ext(combined_ref, dt_used, meta_ref_used, or_used.Q44, no_conj=True))
    else:
        for idx, (dd, rr, oo, dt, mm, mr, ff) in enumerate(zip(data_list,
                                                               ref_list,
                                                               orientation_list,
                                                               dwelltime_list,
                                                               meta_list,
                                                               meta_ref_list,
                                                               files_in)):
            # Add original files to nifti meta information.
            mm.set_standard_def('OriginalFile', [str(ff.name), ])
            fnames_out.append(f'{mainStr}_{idx:03}')
            nifti_mrs_out.append(gen_nifti_mrs_hdr_ext(dd, dt, mm, oo.Q44, no_conj=True))

            if rr is not None:
                mr.set_standard_def('OriginalFile', [str(ff.name), ])
                fnames_out.append(f'{mainStr}_ref_{idx:03}')
                nifti_mrs_out.append(gen_nifti_mrs_hdr_ext(rr, dt, mr, oo.Q44, no_conj=True))

    return nifti_mrs_out, fnames_out


def _process_philips_svs(img, verbose):
    """Process svs data using, using the relevant pathway"""
    if _is_new_format(img):
        return _process_philips_svs_new(img, verbose)
    else:
        return _process_philips_svs_old(img, verbose)


def _process_philips_svs_old(img, verbose):
    """Process Philips DICOM SVS data in the 'old' DICOM format"""
    specData = np.frombuffer(img.dcm_data[('2005', '1270')].value, dtype=np.single)
    specDataCmplx = specData[0::2] + 1j * specData[1::2]

    # 1) Extract dicom parameters
    currNiftiOrientation = _enhanced_dcm_svs_to_orientation(img, verbose)

    dwelltime = 1.0 / img.dcm_data[('2005', '1357')].value
    meta = _extractDicomMetadata_old(img)

    return specDataCmplx, None, currNiftiOrientation, dwelltime, meta, None


def _process_philips_svs_new(img, verbose):
    """Process Philips DICOM SVS data in the 'new' DICOM format"""
    specData = np.frombuffer(img.dcm_data[('5600', '0020')].value, dtype=np.single)
    specDataCmplx = specData[0::2] + 1j * specData[1::2]

    # In the one piece of data I have been provided the data is twice as long as indicated (1 avg)
    # and the second half is a water reference.
    # (0028,0008) IS [32]                                               # 2,1 Number of Frames
    nframes = int(img.dcm_data[0x0028, 0x0008].value)
    spec_points = img.dcm_data.SpectroscopyAcquisitionDataColumns
    reference_inlcuded = img.dcm_data[0x0018, 0x9199].value == 'YES'

    npoints_expected = spec_points * nframes
    if npoints_expected != specDataCmplx.size:
        raise IncorrectDataDimensionsError(
            f'Number of points {specDataCmplx.size}, does not equal frames ({nframes}) '
            f'* spectral points ({spec_points}).')

    if reference_inlcuded:
        spec_data_main = specDataCmplx[:spec_points]
        spec_data_ref = specDataCmplx[spec_points:]
        spec_data_main = spec_data_main.reshape(nframes - 1, spec_points).T
        spec_data_main = spec_data_main.squeeze()
    else:
        spec_data_main = specDataCmplx.reshape(nframes, spec_points).T
        spec_data_main = spec_data_main.squeeze()
        spec_data_ref = None

    # 1) Extract dicom parameters
    currNiftiOrientation = _enhanced_dcm_svs_to_orientation(img, verbose)

    dwelltime = 1.0 / img.dcm_data.SpectralWidth
    meta = _extractDicomMetadata_new(img)
    if reference_inlcuded:
        meta_r = _extractDicomMetadata_new(img, water_suppressed=False)
    else:
        meta_r = None

    return spec_data_main, spec_data_ref, currNiftiOrientation, dwelltime, meta, meta_r


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
    meta = _extractDicomMetadata_new(img)
    meta_r = _extractDicomMetadata_new(img, water_suppressed=False)

    return spec_data_main, spec_data_ref, currNiftiOrientation, dwelltime, meta, meta_r


def _enhanced_dcm_svs_to_orientation(img, verbose=False):
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
    if _is_new_format(img):
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
    else:
        import warnings
        warnings.warn('The orientation and position code for old format Philips DICOM is untested.\
                       Please provide example data to improve the handling of this format.')
        # As implemented in vespa (vespa/analysis/fileio/dicom_philips.py#L181)
        try:
            section  = img.dcm_data[0x2005, 0x1085][0]
            angle_lr = float(section[0x02005, 0x1056].value)
            angle_ap = float(section[0x02005, 0x1054].value)
            angle_hf = float(section[0x02005, 0x1055].value)
            dim_lr   = float(section[0x02005, 0x1059].value)
            dim_ap   = float(section[0x02005, 0x1057].value)
            dim_hf   = float(section[0x02005, 0x1058].value)
            shift_lr = float(section[0x02005, 0x105c].value)
            shift_ap = float(section[0x02005, 0x105a].value)
            shift_hf = float(section[0x02005, 0x105b].value)

            xyzMM = np.array([dim_lr, dim_ap, dim_hf])

            from scipy.spatial.transform import Rotation
            rot = Rotation.from_euler('xyz', [-angle_lr, -angle_ap, angle_hf], degrees=True)
            # THIS NEEDS TESTING!!! CODE WILL PASS BUT I DON"T TRUST IT AT ALL.
            imageOrientationPatient = rot.as_matrix()[:, :2].T
            imagePositionPatient = [-shift_lr, -shift_ap, shift_hf]

        except KeyError:
            # Default orientation
            print('No orientation infroamtion found in header. Tag (2005, 1085) missing. '
                  'Default orientation will be used.')
            imageOrientationPatient = [[1, 0, 0], [0, 1, 0]]
            imagePositionPatient = [0, 0, 0]
            xyzMM = [float(x) for x in img.dcm_data[0x0028, 0x0030].value] \
                + [float(img.dcm_data[0x0018, 0x0050].value), ]

    currNiftiOrientation = dcm_to_nifti_orientation(imageOrientationPatient,
                                                    imagePositionPatient,
                                                    xyzMM,
                                                    (1, 1, 1),
                                                    half_shift=False,
                                                    verbose=verbose)

    return currNiftiOrientation


def _extractDicomMetadata_old(dcmdata, water_suppressed=True):
    """ Extract information from the nibabel DICOM object to insert into the json header ext.

    Args:
        dcmdata: nibabel.nicom image object
    Returns:
        obj (hdr_ext): NIfTI MRS hdr ext object.
    """

    # Extract required metadata and create hdr_ext object
    try:
        frequency = dcmdata.dcm_data.TransmitterFrequency
    except AttributeError:
        frequency = float(dcmdata.dcm_data[('2001', '1083')].value)

    try:
        nucleus = dcmdata.dcm_data.ResonantNucleus
    except AttributeError:
        nucleus = dcmdata.dcm_data[('2001', '1087')].value

    obj = Hdr_Ext(
        frequency,
        nucleus)

    # private_mrs_tags = dcmdata[('2005', '140f')][0]

    def set_if_present(tag, val):
        if val:
            obj.set_standard_def(tag, val)

    # Standard metadata
    # # 5.1 MRS specific Tags
    # 'EchoTime'
    try:
        echo_time = float(dcmdata.dcm_data.PerFrameFunctionalGroupsSequence[0]
                          .MREchoSequence[0].EffectiveEchoTime) * 1E-3
    except AttributeError:
        echo_time = float(dcmdata.dcm_data[('2005', '1310')].value)
    obj.set_standard_def('EchoTime', echo_time)

    # 'RepetitionTime'
    try:
        rep_time = float(dcmdata.dcm_data.PerFrameFunctionalGroupsSequence[0]
                         .MRTimingAndRelatedParametersSequence[0].RepetitionTime) / 1E3
    except AttributeError:
        rep_time = dcmdata.dcm_data.RepetitionTime
    obj.set_standard_def('RepetitionTime', rep_time)

    # 'InversionTime'
    # Not known
    # 'MixingTime'
    # Not known
    # 'ExcitationFlipAngle'
    try:
        fa = float(dcmdata.dcm_data.PerFrameFunctionalGroupsSequence[0]
                   .MRTimingAndRelatedParametersSequence[0].RepetitionTime)
    except AttributeError:
        fa = dcmdata.dcm_data.FlipAngle
    obj.set_standard_def('ExcitationFlipAngle', fa)

    # 'TxOffset'
    # Not known
    # 'VOI'
    # Not known
    # 'WaterSuppressed'
    obj.set_standard_def('WaterSuppressed', water_suppressed)
    # 'WaterSuppressionType'
    # Not known
    # 'SequenceTriggered'
    # Not known
    # # 5.2 Scanner information
    # 'Manufacturer'
    obj.set_standard_def('Manufacturer', dcmdata.dcm_data.Manufacturer)
    # 'ManufacturersModelName'
    obj.set_standard_def('ManufacturersModelName', dcmdata.dcm_data.ManufacturerModelName)
    # 'DeviceSerialNumber'
    if 'DeviceSerialNumber' in dcmdata.dcm_data:
        obj.set_standard_def('DeviceSerialNumber', str(dcmdata.dcm_data.DeviceSerialNumber))
    # 'SoftwareVersions'
    obj.set_standard_def('SoftwareVersions', str(dcmdata.dcm_data.SoftwareVersions))
    # 'InstitutionName'
    obj.set_standard_def('InstitutionName', dcmdata.dcm_data.InstitutionName)
    # 'InstitutionAddress'
    if 'InstitutionAddress' in dcmdata.dcm_data:
        obj.set_standard_def('InstitutionAddress', dcmdata.dcm_data.InstitutionAddress)
    # 'TxCoil'
    # Not known
    # 'RxCoil'
    # ToDo
    # img.dcm_data.SharedFunctionalGroupsSequence[0].MRReceiveCoilSequence[0].ReceiveCoilName

    # # 5.3 Sequence information
    # 'SequenceName'
    if 'SequenceName' in dcmdata.dcm_data:
        obj.set_standard_def('SequenceName', dcmdata.dcm_data.PulseSequenceName)
    # 'ProtocolName'
    obj.set_standard_def('ProtocolName', dcmdata.dcm_data.ProtocolName)
    # # 5.4 Sequence information
    # 'PatientPosition'
    obj.set_standard_def('PatientPosition', dcmdata.dcm_data.PatientPosition)
    # 'PatientName'
    if dcmdata.dcm_data.PatientName:
        obj.set_standard_def('PatientName', dcmdata.dcm_data.PatientName.family_name)
    # 'PatientID'
    # Not known
    # 'PatientWeight'
    if 'PatientWeight' in dcmdata.dcm_data and dcmdata.dcm_data.PatientWeight:
        obj.set_standard_def('PatientWeight', float(dcmdata.dcm_data.PatientWeight))
    # 'PatientDoB'
    set_if_present('PatientDoB', dcmdata.dcm_data.PatientBirthDate)
    # 'PatientSex'
    set_if_present('PatientSex', dcmdata.dcm_data.PatientSex)

    # # 5.5 Provenance and conversion metadata
    # 'ConversionMethod'
    obj.set_standard_def('ConversionMethod', f'spec2nii v{spec2nii_ver}')
    # 'ConversionTime'
    conversion_time = datetime.now().isoformat(sep='T', timespec='milliseconds')
    obj.set_standard_def('ConversionTime', conversion_time)
    # 'OriginalFile'
    # Added later
    # # 5.6 Spatial information
    # 'kSpace'
    obj.set_standard_def('kSpace', [False, False, False])

    return obj


def _extractDicomMetadata_new(dcmdata, water_suppressed=True):
    """ Extract information from the nibabel DICOM object to insert into the json header ext.

    Args:
        dcmdata: nibabel.nicom image object
    Returns:
        obj (hdr_ext): NIfTI MRS hdr ext object.
    """

    # Extract required metadata and create hdr_ext object
    obj = Hdr_Ext(
        dcmdata.dcm_data.TransmitterFrequency,
        dcmdata.dcm_data.ResonantNucleus)

    def set_if_present(tag, val):
        if val:
            obj.set_standard_def(tag, val)

    # Standard metadata
    # # 5.1 MRS specific Tags
    # 'EchoTime'
    echo_time = float(dcmdata.dcm_data.PerFrameFunctionalGroupsSequence[0]
                      .MREchoSequence[0].EffectiveEchoTime) * 1E-3
    obj.set_standard_def('EchoTime', echo_time)
    # 'RepetitionTime'
    try:
        rep_tim = float(dcmdata.dcm_data.PerFrameFunctionalGroupsSequence[0]
                        .MRTimingAndRelatedParametersSequence[0].RepetitionTime) / 1E3
    except AttributeError:
        rep_tim = float(dcmdata.dcm_data.SharedFunctionalGroupsSequence[0]
                        .MRTimingAndRelatedParametersSequence[0].RepetitionTime) / 1E3

    obj.set_standard_def('RepetitionTime', rep_tim)
    # 'InversionTime'
    # Not known
    # 'MixingTime'
    # Not known
    # 'ExcitationFlipAngle'
    try:
        fa = float(dcmdata.dcm_data.PerFrameFunctionalGroupsSequence[0]
                   .MRTimingAndRelatedParametersSequence[0].FlipAngle)
    except AttributeError:
        fa = float(dcmdata.dcm_data.SharedFunctionalGroupsSequence[0]
                   .MRTimingAndRelatedParametersSequence[0].FlipAngle)
    obj.set_standard_def('ExcitationFlipAngle', fa)
    # 'TxOffset'
    # Not known
    # 'VOI'
    # Not known
    # 'WaterSuppressed'
    obj.set_standard_def('WaterSuppressed', water_suppressed)
    # 'WaterSuppressionType'
    # Not known
    # 'SequenceTriggered'
    # Not known
    # # 5.2 Scanner information
    # 'Manufacturer'
    obj.set_standard_def('Manufacturer', dcmdata.dcm_data.Manufacturer)
    # 'ManufacturersModelName'
    obj.set_standard_def('ManufacturersModelName', dcmdata.dcm_data.ManufacturerModelName)
    # 'DeviceSerialNumber'
    if 'DeviceSerialNumber' in dcmdata.dcm_data:
        obj.set_standard_def('DeviceSerialNumber', str(dcmdata.dcm_data.DeviceSerialNumber))
    # 'SoftwareVersions'
    obj.set_standard_def('SoftwareVersions', str(dcmdata.dcm_data.SoftwareVersions))
    # 'InstitutionName'
    if 'InstitutionName' in dcmdata.dcm_data:
        obj.set_standard_def('InstitutionName', dcmdata.dcm_data.InstitutionName)
    # 'InstitutionAddress'
    if 'InstitutionAddress' in dcmdata.dcm_data:
        obj.set_standard_def('InstitutionAddress', dcmdata.dcm_data.InstitutionAddress)
    # 'TxCoil'
    # Not known
    # 'RxCoil'
    # ToDo
    # img.dcm_data.SharedFunctionalGroupsSequence[0].MRReceiveCoilSequence[0].ReceiveCoilName

    # # 5.3 Sequence information
    # 'SequenceName'
    obj.set_standard_def('SequenceName', dcmdata.dcm_data.PulseSequenceName)
    # 'ProtocolName'
    if 'ProtocolName' in dcmdata.dcm_data:
        obj.set_standard_def('ProtocolName', dcmdata.dcm_data.ProtocolName)
    # # 5.4 Sequence information
    # 'PatientPosition'
    obj.set_standard_def('PatientPosition', dcmdata.dcm_data.PatientPosition)
    # 'PatientName'
    if dcmdata.dcm_data.PatientName:
        obj.set_standard_def('PatientName', dcmdata.dcm_data.PatientName.family_name)
    # 'PatientID'
    # Not known
    # 'PatientWeight'
    if 'PatientWeight' in dcmdata.dcm_data and dcmdata.dcm_data.PatientWeight:
        obj.set_standard_def('PatientWeight', float(dcmdata.dcm_data.PatientWeight))
    # 'PatientDoB'
    set_if_present('PatientDoB', dcmdata.dcm_data.PatientBirthDate)
    # 'PatientSex'
    set_if_present('PatientSex', dcmdata.dcm_data.PatientSex)

    # # 5.5 Provenance and conversion metadata
    # 'ConversionMethod'
    obj.set_standard_def('ConversionMethod', f'spec2nii v{spec2nii_ver}')
    # 'ConversionTime'
    conversion_time = datetime.now().isoformat(sep='T', timespec='milliseconds')
    obj.set_standard_def('ConversionTime', conversion_time)
    # 'OriginalFile'
    # Added later
    # # 5.6 Spatial information
    # 'kSpace'
    obj.set_standard_def('kSpace', [False, False, False])

    return obj
