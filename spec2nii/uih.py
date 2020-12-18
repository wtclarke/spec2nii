"""spec2nii module containing functions specific to interpreting UIH DICOM
Possibly also suitable for other manufacturer's DICOM (except Siemens)
Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford
"""
import numpy as np
import nibabel.nicom.dicomwrappers
from spec2nii.dcm2niiOrientation.orientationFuncs import nifti_dicom2mat
from spec2nii.nifti_orientation import NIFTIOrient
from spec2nii import nifti_mrs
from datetime import datetime


def svs_or_CSI(img):
    """Identify from the headers whether data is CSI or SVS."""
    # Currently this only looks in-plane. I don't have examples of CSI from
    # UIH.
    if np.product(img.image_shape) > 1.0:
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
    mainStr = ''
    for idx, fn in enumerate(files_in):
        if verbose:
            print(f'Converting dicom file {fn}')

        img = nibabel.nicom.dicomwrappers.wrapper_from_file(fn)

        mrs_type = svs_or_CSI(img)

        if mrs_type == 'SVS':
            specDataCmplx, orientation, dwelltime, meta_obj = process_uih_svs(img)

            newshape = (1, 1, 1) + specDataCmplx.shape
            specDataCmplx = specDataCmplx.reshape(newshape)

        else:
            specDataCmplx, orientation, dwelltime, meta_obj = process_uih_csi(img)

        data_list.append(specDataCmplx)
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

        dt_used = dwelltime_list[0]
        or_used = orientation_list[0]

        # Add original files to nifti meta information.
        meta_used = meta_list[0]
        meta_used.set_standard_def('OriginalFile', [str(ff) for ff in files_in])

        # Combine data into 5th dimension if needed
        combined_data = np.stack(data_list, axis=-1)

        # Add dimension information (if not None for default)
        if tag:
            meta_used.set_dim_info(0, tag)

        # Create NIFTI MRS object.
        nifti_mrs_out.append(nifti_mrs.NIfTI_MRS(combined_data, or_used.Q44, dt_used, meta_used))
    else:
        for idx, (dd, oo, dt, mm, ff) in enumerate(zip(data_list,
                                                   orientation_list,
                                                   dwelltime_list,
                                                   meta_list,
                                                   files_in)):
            # Add original files to nifti meta information.
            mm.set_standard_def('OriginalFile', [str(ff), ])
            fnames_out.append(f'{mainStr}_{idx:03}')
            nifti_mrs_out.append(nifti_mrs.NIfTI_MRS(dd, oo.Q44, dt, mm))

    return nifti_mrs_out, fnames_out


def process_uih_svs(img):
    """Process UIH DICOM SVS data"""

    specData = np.frombuffer(img.dcm_data[('5600', '0020')].value, dtype=np.single)
    specDataCmplx = specData[0::2] - 1j * specData[1::2]

    # 1) Extract dicom parameters
    imageOrientationPatient = img.image_orient_patient.T
    imagePositionPatient = img.image_position
    xyzMM = np.asarray(img.voxel_sizes)

    # 2) in style of dcm2niix
    # a) calculate Q44
    Q44 = nifti_dicom2mat(imageOrientationPatient, imagePositionPatient, xyzMM)
    # b) calculate nifti quaternion parameters
    Q44[:2, :] *= -1
    # 3) place in data class for nifti orientation parameters
    currNiftiOrientation = NIFTIOrient(Q44)
    dwelltime = 1.0 / img.dcm_data.SpectralWidth
    meta = extractDicomMetadata(img)

    return specDataCmplx, currNiftiOrientation, dwelltime, meta


def process_uih_csi(img):
    raise NotImplementedError('UIH CSI not currently handled. No test data availible!')


def extractDicomMetadata(dcmdata):
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
    obj.set_standard_def('DeviceSerialNumber', str(dcmdata.dcm_data.DeviceSerialNumber))
    obj.set_standard_def('SoftwareVersions', dcmdata.dcm_data.SoftwareVersions)

    obj.set_standard_def('InstitutionName', dcmdata.dcm_data.InstitutionName)
    obj.set_standard_def('InstitutionAddress', dcmdata.dcm_data.InstitutionAddress)

    # No apparent Rx Coil name in header.

    # Some sequence information
    obj.set_standard_def('SequenceName', dcmdata.dcm_data.SequenceName)
    obj.set_standard_def('ProtocolName', dcmdata.dcm_data.ProtocolName)

    # Some subject information
    obj.set_standard_def('PatientPosition', dcmdata.dcm_data.PatientPosition)

    def set_if_present(tag, val):
        if val:
            obj.set_standard_def(tag, val)
    set_if_present('PatientDoB', dcmdata.dcm_data.PatientBirthDate)
    set_if_present('PatientSex', dcmdata.dcm_data.PatientSex)

    if dcmdata.dcm_data.PatientName:
        obj.set_standard_def('PatientName', dcmdata.dcm_data.PatientName.family_name)

    if dcmdata.dcm_data.PatientWeight:
        obj.set_standard_def('PatientWeight', float(dcmdata.dcm_data.PatientWeight))

    # Timing and sequence parameters
    obj.set_standard_def('EchoTime', float(dcmdata.dcm_data.EchoTime) * 1E-3)
    if 'InversionTime' in dcmdata.dcm_data:
        obj.set_standard_def('InversionTime', float(dcmdata.dcm_data.InversionTime) * 1E-3)
    obj.set_standard_def('ExcitationFlipAngle', float(dcmdata.dcm_data.FlipAngle))
    obj.set_standard_def('RepetitionTime', float(dcmdata.dcm_data.RepetitionTime) / 1E3)
    # TO DO  - nibabel might need updating.
    # obj.set_standard_def('TxOffset', )

    # Conversion information
    obj.set_standard_def('ConversionMethod', 'spec2nii')
    conversion_time = datetime.now().isoformat(sep='T', timespec='milliseconds')
    obj.set_standard_def('ConversionTime', conversion_time)

    return obj
