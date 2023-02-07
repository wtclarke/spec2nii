"""spec2nii module containing functions specific to interpreting Siemens DICOM
Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford
"""
from datetime import datetime
import re
import warnings

import numpy as np
import nibabel.nicom.dicomwrappers
from nibabel.nicom import csareader as csar

from mapvbvd.read_twix_hdr import parse_buffer

from nifti_mrs.create_nmrs import gen_nifti_mrs_hdr_ext
from nifti_mrs.hdr_ext import Hdr_Ext

from spec2nii.dcm2niiOrientation.orientationFuncs import dcm_to_nifti_orientation
from spec2nii import __version__ as spec2nii_ver


class inconsistentDataError(Exception):
    pass


class missingTagError(Exception):
    pass


class IncompatibleSOPClassUID(Exception):
    pass


def xa_or_vx(img):
    """Determine siemens DICOM type

    :param img: DICOM image object
    :type img: nibabel.nicom.dicomwrappers.SiemensWrapper
    :raises IncompatibleSOPClassUID: Raised if not a Siemens Syngo Non Image Storage or MRSpectroscopyStorage dicom file
    :return: String either 'xa' or 'vx'
    :rtype: str
    """
    if img.dcm_data.SOPClassUID == '1.3.12.2.1107.5.9.1':
        return 'vx'
    elif img.dcm_data.SOPClassUID == '1.2.840.10008.5.1.4.1.1.4.2':
        return 'xa'
    else:
        if img.dcm_data.SOPClassUID == '1.2.840.10008.5.1.4.1.1.4':
            raise IncompatibleSOPClassUID(
                f'spec2nii detected SOPClassUID {img.dcm_data.SOPClassUID}.'
                ' This normaly contains MR imaging (not spectroscopy) data.'
                ' This data was collected on a'
                f' {img.dcm_data.SoftwareVersions} baseline scanner.'
                ' spec2nii is tested on VA-VE, XA20, and XA30 DICOM files.')
        else:
            raise IncompatibleSOPClassUID(
                'spec2nii does not recognise this SOPClassUID '
                f'{img.dcm_data.SOPClassUID} as MRS data. This data was collected on a'
                f' {img.dcm_data.SoftwareVersions} baseline scanner.'
                ' spec2nii is tested on VA-VE, XA20, and XA30 DICOM files.')


def svs_or_CSI(img):
    """Identify from the csa headers whether data is CSI or SVS."""
    if xa_or_vx(img) == 'vx':
        rows = img.csa_header['tags']['Rows']['items'][0]
        cols = img.csa_header['tags']['Columns']['items'][0]
        slices = img.csa_header['tags']['NumberOfFrames']['items'][0]
    elif xa_or_vx(img) == 'xa':
        spec_fov_geom = img.dcm_data.SharedFunctionalGroupsSequence[0].MRSpectroscopyFOVGeometrySequence[0]
        rows = spec_fov_geom.SpectroscopyAcquisitionPhaseRows
        cols = spec_fov_geom.SpectroscopyAcquisitionPhaseColumns
        slices = spec_fov_geom.SpectroscopyAcquisitionOutOfPlanePhaseSteps

    if np.prod([rows, cols, slices]) > 1.0:
        return 'CSI'
    else:
        return 'SVS'


def multi_file_dicom(files_in, fname_out, tag, verbose):
    """Parse a list of Siemens DICOM files"""

    # Convert each file (combine after)
    data_list = []
    orientation_list = []
    dwelltime_list = []
    meta_list = []
    series_num = []
    inst_num = []
    reference = []
    str_suffix = []
    mainStr = ''
    for idx, fn in enumerate(files_in):
        if verbose:
            print(f'Converting dicom file {fn}')

        img = nibabel.nicom.dicomwrappers.wrapper_from_file(fn)

        try:
            mrs_type = svs_or_CSI(img)
        except IncompatibleSOPClassUID as exc:
            if len(files_in) == 1:
                raise exc
            else:
                print(f'Skipping {fn}.')
                print('Raised IncompatibleSOPClassUID error. Moving to next file.')
                print(f'Message: {str(exc)}\n')
                continue

        if mrs_type == 'SVS':
            specDataCmplx, orientation, dwelltime, meta_obj = process_siemens_svs(img, verbose=verbose)

            newshape = (1, 1, 1) + specDataCmplx.shape
            specDataCmplx = specDataCmplx.reshape(newshape)

        else:
            specDataCmplx, orientation, dwelltime, meta_obj = process_siemens_csi(img, verbose=verbose)

        data_list.append(specDataCmplx)
        orientation_list.append(orientation)
        dwelltime_list.append(dwelltime)
        meta_list.append(meta_obj)

        series_num.append(int(img.dcm_data.SeriesNumber))
        inst_num.append(int(img.dcm_data.InstanceNumber))

        ref_ind, str_suf = identify_integrated_references(img, img.dcm_data.InstanceNumber)
        reference.append(ref_ind)
        str_suffix.append(str_suf)

        if idx == 0:
            if fname_out:
                mainStr = fname_out
            elif 'SeriesDescription' in img.dcm_data:
                mainStr = img.dcm_data.SeriesDescription
            elif 'SeriesInstanceUID' in img.dcm_data:
                mainStr = img.dcm_data.SeriesInstanceUID
            else:
                raise missingTagError("Neither SeriesDescription or SeriesInstanceUID tags defined."
                                      " Please specify an output filename using '-f'")

    # Sort by series and instance number
    data_list = np.asarray(data_list)
    orientation_list = np.asarray(orientation_list)
    dwelltime_list = np.asarray(dwelltime_list)
    meta_list = np.asarray(meta_list)
    series_num = np.asarray(series_num)
    inst_num = np.asarray(inst_num)
    reference = np.asarray(reference)
    str_suffix = np.asarray(str_suffix)
    files_in = np.asarray(files_in)

    sort_index = np.lexsort((inst_num, series_num))  # Sort by series then inst

    data_list = data_list[sort_index, :]
    orientation_list = orientation_list[sort_index]
    dwelltime_list = dwelltime_list[sort_index]
    meta_list = meta_list[sort_index]
    series_num = series_num[sort_index]
    inst_num = inst_num[sort_index]
    reference = reference[sort_index]
    str_suffix = str_suffix[sort_index]
    files_in = files_in[sort_index]

    group_ind = []
    for sn in np.unique(series_num):
        for rn in np.unique(reference):
            group_ind.append(list(np.where(np.logical_and(series_num == sn, reference == rn))[0]))

    if verbose:
        print(f'Sorted series numbers: {series_num}')
        print(f'Sorted instance numbers: {inst_num}')
        print(f'Sorted reference index: {reference}')
        print(f'Output groups: {group_ind}')

    nifti_mrs_out, fnames_out = [], []
    for idx, gr in enumerate(group_ind):

        # If data shape, orientation, dwelltime match then
        # proceed
        def not_equal(lst):
            return lst[:-1] != lst[1:]

        if not_equal([d.shape for d in data_list[gr]])\
                and not_equal([o.Q44.tolist() for o in orientation_list[gr]])\
                and not_equal(dwelltime_list[gr]):
            raise inconsistentDataError('Shape, orientation and dwelltime must match in combined data.')

        fnames_out.append(mainStr + str_suffix[gr[0]])

        dt_used = dwelltime_list[gr[0]]
        or_used = orientation_list[gr[0]]

        # Add original files to nifti meta information.
        meta_used = meta_list[gr[0]]
        meta_used.set_standard_def('OriginalFile', [str(ff) for ff in files_in[gr]])

        # Combine data into 5th dimension if needed
        data_in_gr = data_list[gr]
        if len(data_in_gr) > 1:
            combined_data = np.stack(data_in_gr, axis=-1)
        else:
            combined_data = data_in_gr[0]

        # Add dimension information (if not None for default)
        if tag:
            meta_used.set_dim_info(0, tag)

        # Create NIFTI MRS object.
        try:
            nifti_mrs_out.append(gen_nifti_mrs_hdr_ext(combined_data, dt_used, meta_used, or_used.Q44, no_conj=True))
        except np.linalg.LinAlgError:
            warnings.warn("""The quaternion passed to NIfTI_MRS was singular.
                           Most likely your slice position is not well defined. I have set it to default.""")
            nifti_mrs_out.append(gen_nifti_mrs_hdr_ext(combined_data, dt_used, meta_used, no_conj=True))

    # If there are any identical names then append an index
    seen = np.unique(fnames_out)
    if seen.size < len(fnames_out):
        seen_count = np.zeros(seen.shape, dtype=int)
        fnames_out_checked = []
        for fn in fnames_out:
            if fn in seen:
                seen_index = seen == fn
                fnames_out_checked.append(fn + f'_{seen_count[seen_index][0]:03}')
                seen_count[seen_index] += 1
            else:
                fnames_out_checked.append(fn)

        return nifti_mrs_out, fnames_out_checked
    else:
        return nifti_mrs_out, fnames_out


def process_siemens_svs(img, verbose):
    """Pass to appropriate function for software version."""
    if xa_or_vx(img) == 'vx':
        return process_siemens_svs_vx(img, verbose)
    elif xa_or_vx(img) == 'xa':
        return process_siemens_svs_xa(img, verbose)


def process_siemens_svs_xa(img, verbose):
    """Process Siemens DICOM SVS data acquired on a NumarisX system."""
    specData = np.frombuffer(img.dcm_data[('5600', '0020')].value, dtype=np.single)
    # It appears the phase convention is reversed in the NumarisX systems.
    specDataCmplx = specData[0::2] - 1j * specData[1::2]

    dcm_hdrs = img.dcm_data.PerFrameFunctionalGroupsSequence[0]
    dcm_hdrs1 = img.dcm_data.SharedFunctionalGroupsSequence[0]
    # 1) Extract dicom parameters
    try:
        imageOrientationPatient = np.array(dcm_hdrs.PlaneOrientationSequence[0].ImageOrientationPatient).reshape(2, 3)
    except AttributeError:
        # For VB19 data formatted as MRSpectroscopyStorage this seems to be different to
        # the XA data I've seen. But this could also be a sequence specific thing!
        imageOrientationPatient = np.array(dcm_hdrs1.PlaneOrientationSequence[0].ImageOrientationPatient).reshape(2, 3)
    # VoiPosition - this does not have the FOV shift that imagePositionPatient has
    imagePositionPatient = dcm_hdrs.PlanePositionSequence[0].ImagePositionPatient

    # The first two elements could easily be the reverse of this.
    xyzMM = np.array([float(dcm_hdrs1.PixelMeasuresSequence[0].PixelSpacing[0]),
                      float(dcm_hdrs1.PixelMeasuresSequence[0].PixelSpacing[1]),
                      float(dcm_hdrs1.PixelMeasuresSequence[0].SliceThickness)])

    currNiftiOrientation = dcm_to_nifti_orientation(imageOrientationPatient,
                                                    imagePositionPatient,
                                                    xyzMM,
                                                    (1, 1, 1),
                                                    verbose=verbose)
    dwelltime = 1 / img.dcm_data.SpectralWidth
    meta = extractDicomMetadata_xa(img)

    return specDataCmplx, currNiftiOrientation, dwelltime, meta


def process_siemens_svs_vx(img, verbose):
    """Process Siemens DICOM SVS data acquired on a Numaris4 Vx system."""

    specData = np.frombuffer(img.dcm_data[('7fe1', '1010')].value, dtype=np.single)
    specDataCmplx = specData[0::2] + 1j * specData[1::2]

    # 1) Extract dicom parameters
    imageOrientationPatient = np.array(img.csa_header['tags']['ImageOrientationPatient']['items']).reshape(2, 3)
    # VoiPosition - this does not have the FOV shift that imagePositionPatient has
    imagePositionPatient = img.csa_header['tags']['VoiPosition']['items']
    xyzMM = np.array([img.csa_header['tags']['VoiPhaseFoV']['items'][0],
                      img.csa_header['tags']['VoiReadoutFoV']['items'][0],
                      img.csa_header['tags']['VoiThickness']['items'][0]])

    currNiftiOrientation = dcm_to_nifti_orientation(imageOrientationPatient,
                                                    imagePositionPatient,
                                                    xyzMM,
                                                    (1, 1, 1),
                                                    verbose=verbose)
    dwelltime = img.csa_header['tags']['RealDwellTime']['items'][0] * 1E-9
    meta = extractDicomMetadata_vx(img)

    return specDataCmplx, currNiftiOrientation, dwelltime, meta


def process_siemens_csi(img, verbose):
    """Pass to appropriate CSI function for software version."""
    if xa_or_vx(img) == 'vx':
        return process_siemens_csi_vx(img, verbose)
    elif xa_or_vx(img) == 'xa':
        return process_siemens_csi_xa(img, verbose)


def process_siemens_csi_xa(img, verbose):
    """Process Siemens DICOM CSI data acquired on a NumarisX system."""
    # NOTE: WTC expect to need to reverse the phase convention as above for svs.
    raise NotImplementedError('Method process_siemens_csi_xa not implemented, example data needed!')


def process_siemens_csi_vx(img, verbose):
    """Process Siemens DICOM CSI data acquired on a Numaris4 Vx system."""
    specData = np.frombuffer(img.dcm_data[('7fe1', '1010')].value, dtype=np.single)
    specDataCmplx = specData[0::2] + 1j * specData[1::2]

    rows = img.csa_header['tags']['Rows']['items'][0]
    cols = img.csa_header['tags']['Columns']['items'][0]
    slices = img.csa_header['tags']['NumberOfFrames']['items'][0]
    spectral_points = img.csa_header['tags']['DataPointColumns']['items'][0]

    # I think the above is sufficient, but might need to access more specific Siemens specific parameters
    # In that case I could use the following use the full MrPhoenixProtocol sSpecPara parameters
    # fullcsa = csar.get_csa_header(img.dcm_data, csa_type='series')
    # xprot = parse_buffer(fullcsa['tags']['MrPhoenixProtocol']['items'][0])

    # rows = int(xprot[('sSpecPara', 'lFinalMatrixSizePhase')])
    # cols = int(xprot[('sSpecPara', 'lFinalMatrixSizeRead')])
    # slices = int(xprot[('sSpecPara', 'lFinalMatrixSizeSlice')])
    # spectral_points = int(xprot[('sSpecPara', 'lVectorSize')])

    specDataCmplx = specDataCmplx.reshape((slices, rows, cols, spectral_points))
    specDataCmplx = np.moveaxis(specDataCmplx, (0, 1, 2), (2, 1, 0))

    # 1) Extract dicom parameters
    imageOrientationPatient = np.array(img.csa_header['tags']['ImageOrientationPatient']['items']).reshape(2, 3)
    imagePositionPatient = np.array(img.csa_header['tags']['ImagePositionPatient']['items'])
    xyzMM = np.array([img.csa_header['tags']['PixelSpacing']['items'][0],
                      img.csa_header['tags']['PixelSpacing']['items'][1],
                      img.csa_header['tags']['SliceThickness']['items'][0]])

    # Note that half_shift = True. For an explanation see spec2nii/notes/siemens.md
    currNiftiOrientation = dcm_to_nifti_orientation(imageOrientationPatient,
                                                    imagePositionPatient,
                                                    xyzMM,
                                                    specDataCmplx.shape[:3],
                                                    half_shift=True,
                                                    verbose=verbose)

    dwelltime = img.csa_header['tags']['RealDwellTime']['items'][0] * 1E-9
    meta = extractDicomMetadata_vx(img)

    # Look for a VOI in the data
    meta = _detect_and_fill_voi(img, meta)

    return specDataCmplx, currNiftiOrientation, dwelltime, meta


def _detect_and_fill_voi(img, meta):
    """Look for a VOI volume in the headers of a CSI scan. If present populate the VOI header field.

    :param img: Nibable DICOM image object
    :type img: SiemensWrapper
    :param meta: Existing NIfTI-MRS meta object
    :type meta: hdr_ext
    :return: Modified meta hdr_ext object
    :rtype: hdr_ext
    """
    imageOrientationPatient = np.array(img.csa_header['tags']['ImageOrientationPatient']['items']).reshape(2, 3)
    # VoiPosition - this does not have the FOV shift that imagePositionPatient has
    imagePositionPatient = img.csa_header['tags']['VoiPosition']['items']
    xyzMM = np.array([img.csa_header['tags']['VoiPhaseFoV']['items'][0],
                      img.csa_header['tags']['VoiReadoutFoV']['items'][0],
                      img.csa_header['tags']['VoiThickness']['items'][0]])

    voiNiftiOrientation = dcm_to_nifti_orientation(imageOrientationPatient,
                                                   imagePositionPatient,
                                                   xyzMM,
                                                   (1, 1, 1))

    meta.set_standard_def('VOI', voiNiftiOrientation.Q44.tolist())

    return meta


def extractDicomMetadata_xa(dcmdata):
    """ Extract information from the nibabel DICOM object to insert into the json header ext.
    For Numarix X systems

    Args:
        dcmdata: nibabel.nicom image object
    Returns:
        obj (hdr_ext): NIfTI MRS hdr ext object.
    """

    dcm_hdrs1 = dcmdata.dcm_data.SharedFunctionalGroupsSequence[0]

    # Extract required metadata and create hdr_ext object
    obj = Hdr_Ext(
        dcmdata.dcm_data.TransmitterFrequency,
        dcmdata.dcm_data.ResonantNucleus)

    # Standard defined metadata
    def set_standard_def(nifti_mrs_key, location, key, cast=None):
        try:
            if cast is not None:
                obj.set_standard_def(nifti_mrs_key, cast(getattr(location, key)))
            else:
                obj.set_standard_def(nifti_mrs_key, getattr(location, key))
        except AttributeError:
            pass

    # # 5.1 MRS specific Tags
    # 'EchoTime'
    obj.set_standard_def('EchoTime', float(dcm_hdrs1.MREchoSequence[0].EffectiveEchoTime * 1E-3))
    # 'RepetitionTime'
    obj.set_standard_def(
        'RepetitionTime',
        float(dcm_hdrs1.MRTimingAndRelatedParametersSequence[0].RepetitionTime) / 1E3)
    # 'InversionTime'
    # 'MixingTime'
    # 'ExcitationFlipAngle'
    obj.set_standard_def(
        'ExcitationFlipAngle',
        float(dcm_hdrs1.MRTimingAndRelatedParametersSequence[0].FlipAngle))
    # 'TxOffset'
    # 'VOI'
    # 'WaterSuppressed'
    # 'WaterSuppressionType'
    # 'SequenceTriggered'

    # # 5.2 Scanner information
    # 'Manufacturer'
    set_standard_def('Manufacturer', dcmdata.dcm_data, 'Manufacturer')
    # 'ManufacturersModelName'
    set_standard_def('ManufacturersModelName', dcmdata.dcm_data, 'ManufacturerModelName')
    # 'DeviceSerialNumber'
    set_standard_def('DeviceSerialNumber', dcmdata.dcm_data, 'DeviceSerialNumber', cast=str)
    # 'SoftwareVersions'
    set_standard_def('SoftwareVersions', dcmdata.dcm_data, 'SoftwareVersions')
    # 'InstitutionName'
    set_standard_def('InstitutionName', dcmdata.dcm_data, 'InstitutionName')
    # 'InstitutionAddress'
    set_standard_def('InstitutionAddress', dcmdata.dcm_data, 'InstitutionAddress')
    # 'TxCoil'
    obj.set_standard_def(
        'TxCoil',
        dcm_hdrs1.MRTransmitCoilSequence[0].TransmitCoilName)
    # 'RxCoil'
    obj.set_standard_def(
        'RxCoil',
        dcm_hdrs1.MRReceiveCoilSequence[0].ReceiveCoilName)

    # # 5.3 Sequence information
    # 'SequenceName'
    set_standard_def('SequenceName', dcmdata.dcm_data, 'PulseSequenceName')
    # 'ProtocolName'
    set_standard_def('ProtocolName', dcmdata.dcm_data, 'ProtocolName')
    # # 5.4 Sequence information
    # 'PatientPosition'
    set_standard_def('PatientPosition', dcmdata.dcm_data, 'PatientPosition')
    # 'PatientName'
    set_standard_def('PatientName', dcmdata.dcm_data.PatientName, 'family_name')
    # 'PatientID'
    set_standard_def('PatientID', dcmdata.dcm_data, 'PatientID')
    # 'PatientWeight'
    set_standard_def('PatientWeight', dcmdata.dcm_data, 'PatientWeight', cast=float)
    # 'PatientDoB'
    set_standard_def('PatientDoB', dcmdata.dcm_data, 'PatientBirthDate')
    # 'PatientSex'
    set_standard_def('PatientSex', dcmdata.dcm_data, 'PatientSex')

    # # 5.5 Provenance and conversion metadata
    obj.set_standard_def('ConversionMethod', f'spec2nii v{spec2nii_ver}')
    # 'ConversionTime'
    conversion_time = datetime.now().isoformat(sep='T', timespec='milliseconds')
    obj.set_standard_def('ConversionTime', conversion_time)

    # 'OriginalFile'
    # Set elsewhere

    # # 5.6 Spatial information
    # 'kSpace'
    obj.set_standard_def('kSpace', [False, False, False])

    return obj


def extractDicomMetadata_vx(dcmdata):
    """ Extract information from the nibabel DICOM object to insert into the json header ext.

    Args:
        dcmdata: nibabel.nicom image object
    Returns:
        obj (hdr_ext): NIfTI MRS hdr ext object.
    """

    # Extract required metadata and create hdr_ext object
    obj = Hdr_Ext(
        dcmdata.csa_header['tags']['ImagingFrequency']['items'][0],
        dcmdata.csa_header['tags']['ImagedNucleus']['items'][0])

    # Standard defined metadata
    def set_standard_def(nifti_mrs_key, location, key, cast=None):
        try:
            if cast is not None:
                obj.set_standard_def(nifti_mrs_key, cast(getattr(location, key)))
            else:
                obj.set_standard_def(nifti_mrs_key, getattr(location, key))
        except AttributeError:
            pass

    # # 5.1 MRS specific Tags
    # 'EchoTime' - requires substantial extraction from the full protocol incase there are multiple
    # sub-values, e.g. summing the three TEs of a sLAASER sequence together.
    fullcsa = csar.get_csa_header(dcmdata.dcm_data, csa_type='series')
    xprot = parse_buffer(fullcsa['tags']['MrPhoenixProtocol']['items'][0])
    te_sum = 0
    for teDx in range(128):
        try:
            te_sum += xprot[('alTE', str(teDx))]
        except KeyError:
            break
    obj.set_standard_def('EchoTime', float(te_sum * 1E-6))
    # Original
    # obj.set_standard_def('EchoTime', float(dcmdata.csa_header['tags']['EchoTime']['items'][0] * 1E-3))

    # 'RepetitionTime'
    obj.set_standard_def('RepetitionTime', float(dcmdata.csa_header['tags']['RepetitionTime']['items'][0] / 1E3))
    # 'InversionTime'
    if dcmdata.csa_header['tags']['InversionTime']['n_items'] > 0:
        obj.set_standard_def('InversionTime', float(dcmdata.csa_header['tags']['InversionTime']['items'][0]))
    # 'MixingTime'
    # 'ExcitationFlipAngle'
    obj.set_standard_def('ExcitationFlipAngle', float(dcmdata.csa_header['tags']['FlipAngle']['items'][0]))
    # 'TxOffset'
    # 'VOI'
    # 'WaterSuppressed'
    # 'WaterSuppressionType'
    # 'SequenceTriggered'
    # # 5.2 Scanner information
    # 'Manufacturer'
    set_standard_def('Manufacturer', dcmdata.dcm_data, 'Manufacturer')
    # 'ManufacturersModelName'
    set_standard_def('ManufacturersModelName', dcmdata.dcm_data, 'ManufacturerModelName')
    # 'DeviceSerialNumber'
    set_standard_def('DeviceSerialNumber', dcmdata.dcm_data, 'DeviceSerialNumber', cast=str)
    # 'SoftwareVersions'
    set_standard_def('SoftwareVersions', dcmdata.dcm_data, 'SoftwareVersions')
    # 'InstitutionName'
    set_standard_def('InstitutionName', dcmdata.dcm_data, 'InstitutionName')
    # 'InstitutionAddress'
    set_standard_def('InstitutionAddress', dcmdata.dcm_data, 'InstitutionAddress')
    # 'TxCoil'
    # 'RxCoil'
    if len(dcmdata.csa_header['tags']['ReceivingCoil']['items']) > 0:
        obj.set_standard_def('RxCoil',
                             dcmdata.csa_header['tags']['ReceivingCoil']['items'][0])
    else:
        obj.set_standard_def('RxCoil',
                             dcmdata.csa_header['tags']['ImaCoilString']['items'][0])
    # # 5.3 Sequence information
    # 'SequenceName'
    obj.set_standard_def('SequenceName', dcmdata.csa_header['tags']['SequenceName']['items'][0])
    # 'ProtocolName'
    set_standard_def('ProtocolName', dcmdata.dcm_data, 'ProtocolName')
    # # 5.4 Sequence information
    # 'PatientPosition'
    set_standard_def('PatientPosition', dcmdata.dcm_data, 'PatientPosition')
    # 'PatientName'
    set_standard_def('PatientName', dcmdata.dcm_data.PatientName, 'family_name')
    # 'PatientID'
    # 'PatientWeight'
    set_standard_def('PatientWeight', dcmdata.dcm_data, 'PatientWeight', cast=float)
    # 'PatientDoB'
    set_standard_def('PatientDoB', dcmdata.dcm_data, 'PatientBirthDate')
    # 'PatientSex'
    set_standard_def('PatientSex', dcmdata.dcm_data, 'PatientSex')
    # # 5.5 Provenance and conversion metadata
    obj.set_standard_def('ConversionMethod', f'spec2nii v{spec2nii_ver}')
    # 'ConversionTime'
    conversion_time = datetime.now().isoformat(sep='T', timespec='milliseconds')
    obj.set_standard_def('ConversionTime', conversion_time)
    # 'OriginalFile'
    # Set elsewhere
    # # 5.6 Spatial information
    # 'kSpace'
    obj.set_standard_def('kSpace', [False, False, False])

    # Some additional sequence information
    obj.set_user_def(key='PulseSequenceFile',
                     value=dcmdata.csa_header['tags']['SequenceName']['items'][0],
                     doc='Sequence binary path.')
    # obj.set_user_def(key='IceProgramFile',
    #                  value=mapVBVDHdr['Meas'][('tICEProgramName')],
    #                  doc='Reconstruction binary path.')

    return obj


def identify_integrated_references(img, inst_num):
    '''Heuristics for identifying integrated reference scans in known sequences.
    Sequences handled: CMRR svs_slaserVOI_dkd(2)

    :param img: nibable dicom image

    :return: Ref scan index 0 = not reference scan, higher integer splits into groups.
    :return: name suffix
    '''

    if xa_or_vx(img) == 'xa':
        return 0, ''

    fullcsa = csar.get_csa_header(img.dcm_data, csa_type='series')
    xprot = parse_buffer(fullcsa['tags']['MrPhoenixProtocol']['items'][0])

    # Handle CMRR DKD sequence
    # https://www.cmrr.umn.edu/spectro/
    # SEMI-LASER (MRM 2011, NMB 2019) Release 2016-12
    seq_file_name = xprot[('tSequenceFileName',)].strip('"').lower()
    match = re.search(r'svs_slaservoi_dkd', seq_file_name)
    if match and xprot[('sSpecPara', 'lAutoRefScanMode')] == 8.0:
        num_ref = int(xprot[('sSpecPara', 'lAutoRefScanNo')])
        num_dyn = int(xprot[('lAverages',)])
        total_dyn = num_dyn + (num_ref * 4)
        if inst_num <= num_ref:
            # First ecc calibration references
            return 1, '_ecc'
        elif inst_num <= (num_ref * 2):
            # First quantitation calibration references
            return 2, '_quant'
        elif (total_dyn - (2 * num_ref)) < inst_num <= (total_dyn - num_ref):
            # Second ecc calibration references
            return 1, '_ecc'
        elif (total_dyn - num_ref) < inst_num <= total_dyn:
            # Second quantitation calibration references
            return 2, '_quant'
        else:
            return 0, ''
    else:
        return 0, ''
