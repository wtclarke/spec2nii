"""spec2nii module containing functions specific to interpreting UIH DICOM
Possibly also suitable for other manufacturer's DICOM (except Siemens)
Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford
"""
from datetime import datetime
from warnings import warn

import numpy as np
import nibabel.nicom.dicomwrappers

from nifti_mrs.create_nmrs import gen_nifti_mrs_hdr_ext
from nifti_mrs.hdr_ext import Hdr_Ext

from spec2nii.dcm2niiOrientation.orientationFuncs import dcm_to_nifti_orientation
from spec2nii.nifti_orientation import NIFTIOrient
from spec2nii import __version__ as spec2nii_ver


def svs_or_CSI(img):
    """Identify from the headers whether data is CSI or SVS."""
    # Currently this only looks in-plane. I don't have examples of CSI from
    # UIH.
    if np.prod(img.image_shape) > 1.0:
        return 'CSI'
    else:
        return 'SVS'


def multi_file_dicom(files_in, fname_out, tag, verbose):
    """Parse a list of UIH DICOM files"""

    # Convert each file (combine after)
    data_list = []
    orientation_list = []
    dwelltime_list = []
    meta_list = []
    name_list = []
    mainStr = ''
    for idx, fn in enumerate(files_in):
        if verbose:
            print(f'Converting dicom file {fn}')

        img = nibabel.nicom.dicomwrappers.wrapper_from_file(fn)

        mrs_type = svs_or_CSI(img)

        if mrs_type == 'SVS':
            specDataCmplx, \
                orientation, \
                dwelltime, \
                meta_obj, \
                name_suffix = process_uih_svs(img, verbose)

            # Add original files to nifti meta information.
            for mm in meta_obj:
                mm.set_standard_def('OriginalFile', [str(fn.name), ])

        else:
            specDataCmplx, orientation, dwelltime, meta_obj = process_uih_csi(img, verbose)
            name_suffix = []
            # Add original files to nifti meta information.
            for mm in meta_obj:
                mm.set_standard_def('OriginalFile', [str(fn.name), ])

        data_list.extend(specDataCmplx)
        orientation_list.extend(orientation)
        dwelltime_list.extend(dwelltime)
        meta_list.extend(meta_obj)
        name_list.extend(name_suffix)

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
        and all_equal(dwelltime_list)\
        and all_equal(name_list)

    nifti_mrs_out, fnames_out = [], []
    if combine:
        # Combine files into single MRS NIfTI
        # Single file name
        fnames_out.append(mainStr)

        dt_used = dwelltime_list[0]
        or_used = orientation_list[0]

        # Add original files to nifti meta information.
        meta_used = meta_list[0]
        meta_used.set_standard_def('OriginalFile', [str(ff.name) for ff in files_in])

        # Combine data into 5th dimension if needed
        if len(data_list) > 1:
            combined_data = np.stack(data_list, axis=-1)
        else:
            combined_data = data_list[0]

        # Add dimension information (if not None for default)
        if tag:
            meta_used.set_dim_info(0, tag)

        # Create NIFTI MRS object.
        nifti_mrs_out.append(
            gen_nifti_mrs_hdr_ext(
                combined_data,
                dt_used,
                meta_used,
                or_used.Q44,
                no_conj=True))
    else:
        for idx, (dd, oo, dt, mm, nn) in enumerate(
            zip(
                data_list,
                orientation_list,
                dwelltime_list,
                meta_list,
                name_list)):
            if len(files_in) > 1:
                fnames_out.append(f'{mainStr}_{idx:03}' + nn)
            else:
                fnames_out.append(mainStr + nn)
            nifti_mrs_out.append(
                gen_nifti_mrs_hdr_ext(
                    dd,
                    dt,
                    mm,
                    oo.Q44,
                    no_conj=True))

    return nifti_mrs_out, fnames_out


def process_uih_svs(img, verbose):
    """Process UIH DICOM SVS data"""

    specDataCmplx = []
    meta = []
    dwelltime = []
    orientation = []
    name_suffix = []

    def read_data(tag):
        specData = np.frombuffer(
            img.dcm_data[tag].value,
            dtype=np.single)
        return specData[0::2] + 1j * specData[1::2]

    specDataCmplx.append(read_data(('5600', '0020')))

    # Orientation
    # 1) Extract dicom parameters
    imageOrientationPatient = img.image_orient_patient.T
    imagePositionPatient = img.image_position
    xyzMM = np.asarray(img.voxel_sizes)

    # 2) To nifti
    currNiftiOrientation = dcm_to_nifti_orientation(imageOrientationPatient,
                                                    imagePositionPatient,
                                                    xyzMM,
                                                    (1, 1, 1),
                                                    half_shift=True,
                                                    verbose=verbose)

    base_meta = extractDicomMetadata(img)

    wsupp_meta = base_meta.copy()
    wsupp_meta.set_standard_def('WaterSuppressed', True)
    orientation.append(currNiftiOrientation)
    dwelltime.append(1.0 / img.dcm_data.SpectralWidth)
    meta.append(base_meta.copy())
    name_suffix.append('')

    specDataCmplx.append(read_data(('0065', 'ff07')))
    orientation.append(currNiftiOrientation)
    dwelltime.append(1.0 / img.dcm_data.SpectralWidth)
    wref_meta = base_meta.copy()
    wref_meta.set_standard_def('WaterSuppressed', False)
    meta.append(wref_meta)
    name_suffix.append('_wref')

    if img.dcm_data[('0018', '0024')].value == 'svs_mega':
        name_suffix[0] = '_diff'

        off = read_data(('0065', 'ff11'))
        on = read_data(('0065', 'ff10'))
        specDataCmplx.append(np.stack([off, on], axis=-1))
        name_suffix.append('_off_on')
        orientation.append(currNiftiOrientation)
        dwelltime.append(1.0 / img.dcm_data.SpectralWidth)
        off_on_meta = base_meta.copy()
        off_on_meta.set_standard_def('WaterSuppressed', True)
        off_on_meta.set_dim_info(0, 'DIM_EDIT', hdr={"EditCondition": ["OFF", "ON"]})
        meta.append(off_on_meta)

    specDataCmplx = [data.reshape((1, 1, 1) + data.shape) for data in specDataCmplx]

    return specDataCmplx, orientation, dwelltime, meta, name_suffix


def process_uih_csi(img, verbose):
    """Process UIH DICOM MRSI data"""

    specData = np.frombuffer(img.dcm_data[('5600', '0020')].value, dtype=np.single)
    specDataCmplx = specData[0::2] + 1j * specData[1::2]

    warn('The orientation of UIH MRSI is relatively untested.'
         ' Please contribute data to help fix this!')
    shape = (img.dcm_data.Columns,
             img.dcm_data.Rows,
             int(img.dcm_data.NumberOfFrames),
             img.dcm_data.SpectroscopyAcquisitionDataColumns)

    specDataCmplx = specDataCmplx.reshape(shape)

    # WTC 05/01/21 From the limited test data I have I think there is this transposition.
    specDataCmplx = np.swapaxes(specDataCmplx, 0, 1)

    # Extract dicom parameters
    imageOrientationPatient = img.image_orient_patient.T
    imagePositionPatient = img.image_position
    xyzMM = np.asarray(img.voxel_sizes)

    currNiftiOrientation = dcm_to_nifti_orientation(imageOrientationPatient,
                                                    imagePositionPatient,
                                                    xyzMM,
                                                    specDataCmplx.shape[:3],
                                                    half_shift=True,
                                                    verbose=verbose)

    # UIH specific tweaks
    Q44 = currNiftiOrientation.Q44
    # 1) negate 3rd direction if ('0065', 'ff06') 3rd element is negative
    if float(img.dcm_data[('0065', 'ff06')].value[2]) < 0.0:
        Q44[:3, :3] = Q44[:3, :3] @ np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]])

    # 2) A half voxel shift is needed in the third dimension
    v = np.asarray([0, 0, 0.5]) @ Q44[:3, :3].T
    Q44[:3, 3] += v

    currNiftiOrientation = NIFTIOrient(Q44)

    dwelltime = 1.0 / img.dcm_data.SpectralWidth
    meta = extractDicomMetadata(img)

    return [specDataCmplx], [currNiftiOrientation], [dwelltime], [meta]


def extractDicomMetadata(dcmdata) -> Hdr_Ext:
    """ Extract information from the nibabel DICOM object to insert into the json header ext.

    There seems to be a large 'uprotocol' in tag 0065,1007 but I don't know how to interpret it.

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
    # Standard defined
    # # 5.1 MRS specific Tags
    # 'EchoTime'
    obj.set_standard_def('EchoTime', float(dcmdata.dcm_data.EchoTime) * 1E-3)
    # 'RepetitionTime'
    obj.set_standard_def('RepetitionTime', float(dcmdata.dcm_data.RepetitionTime) / 1E3)
    # 'InversionTime'
    if 'InversionTime' in dcmdata.dcm_data:
        obj.set_standard_def('InversionTime', float(dcmdata.dcm_data.InversionTime) * 1E-3)
    # 'MixingTime'
    # 'ExcitationFlipAngle'
    obj.set_standard_def('ExcitationFlipAngle', float(dcmdata.dcm_data.FlipAngle))
    # 'TxOffset'
    # 'VOI'
    # 'WaterSuppressed'
    # TO DO
    # 'WaterSuppressionType'
    # 'SequenceTriggered'
    # # 5.2 Scanner information
    # 'Manufacturer'
    obj.set_standard_def('Manufacturer', dcmdata.dcm_data.Manufacturer)
    # 'ManufacturersModelName'
    obj.set_standard_def('ManufacturersModelName', dcmdata.dcm_data.ManufacturerModelName)
    # 'DeviceSerialNumber'
    obj.set_standard_def('DeviceSerialNumber', str(dcmdata.dcm_data.DeviceSerialNumber))
    # 'SoftwareVersions'
    obj.set_standard_def('SoftwareVersions', dcmdata.dcm_data.SoftwareVersions)
    # 'InstitutionName'
    obj.set_standard_def('InstitutionName', dcmdata.dcm_data.InstitutionName)
    # 'InstitutionAddress'
    obj.set_standard_def('InstitutionAddress', dcmdata.dcm_data.InstitutionAddress)
    # 'TxCoil'
    # 'RxCoil'
    # No apparent coil name in header.
    # # 5.3 Sequence information
    # 'SequenceName'
    obj.set_standard_def('SequenceName', dcmdata.dcm_data.SequenceName)
    # 'ProtocolName'
    obj.set_standard_def('ProtocolName', dcmdata.dcm_data.ProtocolName)
    # # 5.4 Sequence information
    # 'PatientPosition'
    obj.set_standard_def('PatientPosition', dcmdata.dcm_data.PatientPosition)
    # 'PatientName'
    if dcmdata.dcm_data.PatientName:
        obj.set_standard_def('PatientName', dcmdata.dcm_data.PatientName.family_name)
    # 'PatientID'
    # 'PatientWeight'
    if dcmdata.dcm_data.PatientWeight:
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
    # Set elsewhere
    # # 5.6 Spatial information
    # 'kSpace'
    obj.set_standard_def('kSpace', [False, False, False])

    return obj
