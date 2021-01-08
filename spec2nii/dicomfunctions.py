"""spec2nii module containing functions specific to interpreting Siemens DICOM
Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford
"""
import numpy as np
import nibabel.nicom.dicomwrappers
from spec2nii.dcm2niiOrientation.orientationFuncs import dcm_to_nifti_orientation
from spec2nii import nifti_mrs
from datetime import datetime


def svs_or_CSI(img):
    """Identify from the csa headers whether data is CSI or SVS."""
    rows = img.csa_header['tags']['Rows']['items'][0]
    cols = img.csa_header['tags']['Columns']['items'][0]
    slices = img.csa_header['tags']['NumberOfFrames']['items'][0]

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
    mainStr = ''
    for idx, fn in enumerate(files_in):
        if verbose:
            print(f'Converting dicom file {fn}')

        img = nibabel.nicom.dicomwrappers.wrapper_from_file(fn)

        mrs_type = svs_or_CSI(img)

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

        if idx == 0:
            if fname_out:
                mainStr = fname_out
            else:
                mainStr = img.dcm_data.SeriesDescription

    # If data shape, orientation, dwelltime and series_num match combine
    # into one NIFTI MRS object.
    # Otherwise return a list of files/names
    def all_equal(lst):
        return lst[:-1] == lst[1:]

    combine = all_equal([d.shape for d in data_list])\
        and all_equal([o.Q44.tolist() for o in orientation_list])\
        and all_equal(dwelltime_list)\
        and all_equal(series_num)

    # Sort by series and instance number
    data_list = np.asarray(data_list)
    orientation_list = np.asarray(orientation_list)
    dwelltime_list = np.asarray(dwelltime_list)
    meta_list = np.asarray(meta_list)
    series_num = np.asarray(series_num)
    inst_num = np.asarray(inst_num)

    sort_index = np.lexsort((inst_num, series_num))  # Sort by series then inst

    data_list = data_list[sort_index, :]
    orientation_list = orientation_list[sort_index]
    dwelltime_list = dwelltime_list[sort_index]
    meta_list = meta_list[sort_index]
    series_num = series_num[sort_index]
    inst_num = inst_num[sort_index]

    if verbose:
        print(f'Sorted series numbers: {series_num}')
        print(f'Sorted instance numbers: {inst_num}')

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
        if len(data_list) > 1:
            combined_data = np.stack(data_list, axis=-1)
        else:
            combined_data = data_list[0]

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
            fnames_out.append(f'{mainStr}_{series_num[idx]:03}_{inst_num[idx]:03}')
            nifti_mrs_out.append(nifti_mrs.NIfTI_MRS(dd, oo.Q44, dt, mm))

    return nifti_mrs_out, fnames_out


def process_siemens_svs(img, verbose):
    """Process Siemens DICOM SVS data"""

    specData = np.frombuffer(img.dcm_data[('7fe1', '1010')].value, dtype=np.single)
    specDataCmplx = specData[0::2] - 1j * specData[1::2]

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
    meta = extractDicomMetadata(img)

    return specDataCmplx, currNiftiOrientation, dwelltime, meta


def process_siemens_csi(img, verbose):
    specData = np.frombuffer(img.dcm_data[('7fe1', '1010')].value, dtype=np.single)
    specDataCmplx = specData[0::2] - 1j * specData[1::2]

    rows = img.csa_header['tags']['Rows']['items'][0]
    cols = img.csa_header['tags']['Columns']['items'][0]
    slices = img.csa_header['tags']['NumberOfFrames']['items'][0]
    spectral_points = img.csa_header['tags']['DataPointColumns']['items'][0]

    specDataCmplx = specDataCmplx.reshape((slices, rows, cols, spectral_points))
    specDataCmplx = np.moveaxis(specDataCmplx, (0, 1, 2), (2, 1, 0))

    # 1) Extract dicom parameters
    imageOrientationPatient = np.array(img.csa_header['tags']['ImageOrientationPatient']['items']).reshape(2, 3)
    imagePositionPatient = np.array(img.csa_header['tags']['ImagePositionPatient']['items'])
    xyzMM = np.array([img.csa_header['tags']['PixelSpacing']['items'][0],
                      img.csa_header['tags']['PixelSpacing']['items'][1],
                      img.csa_header['tags']['SliceThickness']['items'][0]])

    # Note that half_shift = True. For an explination see spec2nii/notes/seimens.md
    currNiftiOrientation = dcm_to_nifti_orientation(imageOrientationPatient,
                                                    imagePositionPatient,
                                                    xyzMM,
                                                    specDataCmplx.shape[:3],
                                                    half_shift=True,
                                                    verbose=verbose)

    dwelltime = img.csa_header['tags']['RealDwellTime']['items'][0] * 1E-9
    meta = extractDicomMetadata(img)

    return specDataCmplx, currNiftiOrientation, dwelltime, meta


def extractDicomMetadata(dcmdata):
    """ Extract information from the nibabel DICOM object to insert into the json header ext.

    Args:
        dcmdata: nibabel.nicom image object
    Returns:
        obj (hdr_ext): NIfTI MRS hdr ext object.
    """

    # Extract required metadata and create hdr_ext object
    obj = nifti_mrs.hdr_ext(dcmdata.csa_header['tags']['ImagingFrequency']['items'][0],
                            dcmdata.csa_header['tags']['ImagedNucleus']['items'][0])

    # Some scanner information
    def set_standard_def(nifti_mrs_key, location, key, cast=None):
        try:
            if cast is not None:
                obj.set_standard_def(nifti_mrs_key, cast(getattr(location, key)))
            else:
                obj.set_standard_def(nifti_mrs_key, getattr(location, key))
        except AttributeError:
            pass

    set_standard_def('Manufacturer', dcmdata.dcm_data, 'Manufacturer')
    set_standard_def('ManufacturersModelName', dcmdata.dcm_data, 'ManufacturerModelName')
    set_standard_def('DeviceSerialNumber', dcmdata.dcm_data, 'DeviceSerialNumber', cast=str)
    set_standard_def('SoftwareVersions', dcmdata.dcm_data, 'SoftwareVersions')

    set_standard_def('InstitutionName', dcmdata.dcm_data, 'InstitutionName')
    set_standard_def('InstitutionAddress', dcmdata.dcm_data, 'InstitutionAddress')

    if len(dcmdata.csa_header['tags']['ReceivingCoil']['items']) > 0:
        obj.set_user_def(key='ReceiveCoilName',
                         value=dcmdata.csa_header['tags']['ReceivingCoil']['items'][0],
                         doc='Rx coil name.')
    else:
        obj.set_user_def(key='ReceiveCoilName',
                         value=dcmdata.csa_header['tags']['ImaCoilString']['items'][0],
                         doc='Rx coil name.')

    # Some sequence information
    obj.set_standard_def('SequenceName', dcmdata.csa_header['tags']['SequenceName']['items'][0])
    set_standard_def('ProtocolName', dcmdata.dcm_data, 'ProtocolName')

    obj.set_user_def(key='PulseSequenceFile',
                     value=dcmdata.csa_header['tags']['SequenceName']['items'][0],
                     doc='Sequence binary path.')
    # obj.set_user_def(key='IceProgramFile',
    #                  value=mapVBVDHdr['Meas'][('tICEProgramName')],
    #                  doc='Reconstruction binary path.')

    # Some subject information
    set_standard_def('PatientPosition', dcmdata.dcm_data, 'PatientPosition')
    set_standard_def('PatientName', dcmdata.dcm_data.PatientName, 'family_name')
    set_standard_def('PatientWeight', dcmdata.dcm_data, 'PatientWeight', cast=float)
    set_standard_def('PatientDoB', dcmdata.dcm_data, 'PatientBirthDate')
    set_standard_def('PatientSex', dcmdata.dcm_data, 'PatientSex')

    # Timing and sequence parameters
    obj.set_standard_def('EchoTime', dcmdata.csa_header['tags']['EchoTime']['items'][0] * 1E-3)
    if dcmdata.csa_header['tags']['InversionTime']['n_items'] > 0:
        obj.set_standard_def('InversionTime', dcmdata.csa_header['tags']['InversionTime']['items'][0])
    obj.set_standard_def('ExcitationFlipAngle', dcmdata.csa_header['tags']['FlipAngle']['items'][0])
    obj.set_standard_def('RepetitionTime', dcmdata.csa_header['tags']['RepetitionTime']['items'][0] / 1E3)
    # TO DO  - nibabel might need updating.
    # obj.set_standard_def('TxOffset', )

    # Conversion information
    obj.set_standard_def('ConversionMethod', 'spec2nii')
    conversion_time = datetime.now().isoformat(sep='T', timespec='milliseconds')
    obj.set_standard_def('ConversionTime', conversion_time)

    return obj
