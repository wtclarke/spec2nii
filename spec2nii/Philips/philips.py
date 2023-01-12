"""spec2nii module containing functions specific to interpreting Philips formats
Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford
"""
from datetime import datetime
from ast import literal_eval

import numpy as np
from nifti_mrs.create_nmrs import gen_nifti_mrs_hdr_ext
from nifti_mrs.hdr_ext import Hdr_Ext

from spec2nii.nifti_orientation import NIFTIOrient, calc_affine
from spec2nii import __version__ as spec2nii_ver


def read_sdat_spar_pair(sdat_file, spar_file, shape=None, tags=None, fileout=None, special=None):
    """Read and convert SPAR/SDAT pairs from Philips scanners

    :param sdat_file: Path to .SDAT file
    :type sdat_file: pathlib.Path
    :param spar_file: Path to .SPAR file
    :type spar_file: pathlib.Path
    :param shape: List of dimension shapes to reshape array, defaults to None
    :type shape: list of ints, optional
    :param tags: List of higher dimension tags, defaults to None
    :type tags: List of strings, optional
    :param fileout: File name string passed by user, defaults to None/stem of file
    :type fileout: str, optional
    :param special: Identifier for special-cased sequence, defaults to None
    :type special: str, optional
    :return: List of NIFTI-MRS objects
    :rtype: list of nifti_mrs.nifti_mrs.NIFTI_MRS
    :return: List of file name parts
    :rtype: list of strings
    """
    spar_params = read_spar(spar_file)
    data = read_sdat(sdat_file,
                     spar_params['samples'],
                     spar_params['rows'])

    if data.ndim < 4:
        data = data.reshape((1, 1, 1) + data.shape)

    if shape is not None:
        it_shape = data.shape[:4]
        data = data.reshape(it_shape + tuple(shape), order='F')

    # Move to right handed frame
    data = data.conj()

    # Dwelltime
    dwelltime = 1.0 / float(spar_params["sample_frequency"])

    # Meta
    meta = spar_to_nmrs_hdrext(spar_params)
    meta.set_standard_def('OriginalFile', [sdat_file.name])

    for idx, tag in enumerate(tags):
        if tag is not None:
            meta.set_dim_info(idx, tag)

    # Orientation
    if spar_params["volume_selection_enable"] == "yes":
        affine = _philips_orientation(spar_params)
    else:
        # Use default affine
        affine = np.diag(np.array([10000, 10000, 10000, 1]))
    orientation = NIFTIOrient(affine)

    # name of output
    if fileout is not None:
        mainStr = fileout
    else:
        mainStr = sdat_file.stem

    # Special cases
    if spar_params['scan_id'].lower() == 'hyper'\
            or special.lower() == 'hyper':
        return _special_case_hyper(data, dwelltime, meta, orientation.Q44, mainStr)
    else:
        # Normal case
        return [gen_nifti_mrs_hdr_ext(data, dwelltime, meta, orientation.Q44, no_conj=True), ],\
               [mainStr, ]


def read_spar(filename):
    '''Read the .spar file.
    :param filename: file path

    :return: dict of parameters read from spar file
    :rtype: dict
    '''

    parameter_dict = {}
    with open(filename, 'r') as f:
        for line in f:
            # ignore comments (!) and empty lines
            if line == "\n" or line.startswith("!"):
                continue

            # Handle
            key, value = map(str.strip, line.split(":", 1))
            try:
                val = literal_eval(value)
            except (ValueError, SyntaxError):
                if value == '':
                    val = None
                else:
                    val = value

            parameter_dict.update({key: val})

    return parameter_dict


def read_sdat(filename, samples, rows):
    '''Read the .sdat file.
    :param filename: File path
    :param samples: Number of spectral points
    :param rows: Number of rows of data
    '''
    with open(filename, 'rb') as f:
        raw = f.read()

    floats = _vax_to_ieee_single_float(raw)
    data_iter = iter(floats)
    complex_iter = (complex(r, i) for r, i in zip(data_iter, data_iter))
    raw_data = np.fromiter(complex_iter, "complex64")
    raw_data = np.reshape(raw_data, (rows, samples)).T.squeeze()

    return raw_data


def _philips_orientation(params):
    '''Calculate the orientation affine from the spar parameters.'''

    angle_lr = params['lr_angulation']
    angle_ap = params['ap_angulation']
    angle_hf = params['cc_angulation']
    angles = [-angle_lr, -angle_ap, angle_hf]

    dim_lr = params['lr_size']
    dim_ap = params['ap_size']
    dim_hf = params['cc_size']
    dimensions = [dim_lr, dim_ap, dim_hf]

    shift_lr = params['lr_off_center']
    shift_ap = params['ap_off_center']
    shift_hf = params['cc_off_center']
    shift = [-shift_lr, -shift_ap, shift_hf]

    return calc_affine(angles, dimensions, shift)


def spar_to_nmrs_hdrext(spar_dict):
    """ Extract information from the dict of keys read from the spar file to insert into the json header ext.

    :param dict spar_dict: key-value pairs read from spar file
    :return: NIfTI MRS hdr ext object.
    """

    # Extract required metadata and create hdr_ext object
    cf = float(spar_dict["synthesizer_frequency"]) / 1E6
    obj = Hdr_Ext(
        cf,
        spar_dict["nucleus"])

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
    obj.set_standard_def('EchoTime', float(spar_dict['echo_time']) * 1E-3)
    # 'RepetitionTime'
    obj.set_standard_def('RepetitionTime', float(spar_dict['repetition_time'] / 1E3))
    # 'InversionTime'
    if spar_dict['spectrum_inversion_time'] > 0:
        obj.set_standard_def('InversionTime', float(spar_dict['spectrum_inversion_time']))
    # 'MixingTime'
    # 'ExcitationFlipAngle'
    # 'TxOffset'
    obj.set_standard_def('TxOffset', float(spar_dict['offset_frequency'] / cf))
    # 'VOI'
    # 'WaterSuppressed'
    # No apparent parameter stored in the SPAR info.
    # 'WaterSuppressionType'
    # 'SequenceTriggered'
    # # 5.2 Scanner information
    # 'Manufacturer'
    obj.set_standard_def('Manufacturer', 'Philips')
    # 'ManufacturersModelName'
    # 'DeviceSerialNumber'
    # 'SoftwareVersions'
    set_standard_def('SoftwareVersions', spar_dict, 'equipment_sw_verions')
    # 'InstitutionName'
    # 'InstitutionAddress'
    # 'TxCoil'
    # 'RxCoil'
    # # 5.3 Sequence information
    # 'SequenceName'
    # 'ProtocolName'
    set_standard_def('ProtocolName', spar_dict, 'scan_id')
    # # 5.4 Sequence information
    # 'PatientPosition'
    try:
        obj.set_standard_def('PatientPosition', spar_dict['patient_position'] + ' ' + spar_dict['patient_orientation'])
    except AttributeError:
        pass
    # 'PatientName'
    set_standard_def('PatientName', spar_dict, 'patient_name')
    # 'PatientID'
    # 'PatientWeight'
    # 'PatientDoB'
    set_standard_def('PatientDoB', spar_dict, 'patient_birth_date')
    # 'PatientSex'
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


# From VESPA - BSD license.
def _vax_to_ieee_single_float(data):
    """Converts a float in Vax format to IEEE format.

    data should be a single string of chars that have been read in from
    a binary file. These will be processed 4 at a time into float values.
    Thus the total number of byte/chars in the string should be divisible
    by 4.

    Based on VAX data organization in a byte file, we need to do a bunch of
    bitwise operations to separate out the numbers that correspond to the
    sign, the exponent and the fraction portions of this floating point
    number

    role :      S        EEEEEEEE      FFFFFFF      FFFFFFFF      FFFFFFFF
    bits :      1        2      9      10                               32
    bytes :     byte2           byte1               byte4         byte3

    This is taken from the VESPA project source code under a BSD licence.
    """
    f = []
    nfloat = int(len(data) / 4)
    for i in range(nfloat):

        byte2 = data[0 + i * 4]
        byte1 = data[1 + i * 4]
        byte4 = data[2 + i * 4]
        byte3 = data[3 + i * 4]

        # hex 0x80 = binary mask 10000000
        # hex 0x7f = binary mask 01111111

        sign = (byte1 & 0x80) >> 7
        expon = ((byte1 & 0x7f) << 1) + ((byte2 & 0x80) >> 7)
        fract = ((byte2 & 0x7f) << 16) + (byte3 << 8) + byte4

        if sign == 0:
            sign_mult = 1.0
        else:
            sign_mult = -1.0

        if 0 < expon:
            # note 16777216.0 == 2^24
            val = sign_mult * (0.5 + (fract / 16777216.0)) * pow(2.0, expon - 128.0)
            f.append(val)
        elif expon == sign == 0:
            f.append(0)
        else:
            f.append(0)
            # may want to raise an exception here ...

    return f


def _special_case_hyper(data, dwelltime, meta, orientation, fout_str):
    # Reorganise the data. This unfortunately makes hardcoded assumptions about the size of each part.
    data_short_te = data[:, :, :, :, :32]
    data_edited = data[:, :, :, :, 32:]
    data_edited = data_edited.T.reshape((56, 4, data.shape[3], 1, 1, 1)).T

    meta_short_te = meta.copy()
    meta_edited = meta.copy()

    edit_pulse_1 = 0.0
    edit_pulse_2 = 0.0
    edit_pulse_off = 0.0
    dim_info = "HERCULES j-difference editing, four conditions"
    dim_header = {"EditCondition": ["A", "B", "C", "D"]}
    edit_pulse_val = {
        "A": {"PulseOffset": [edit_pulse_1, edit_pulse_2], "PulseDuration": 0.02},
        "B": {"PulseOffset": [edit_pulse_off, edit_pulse_2], "PulseDuration": 0.02},
        "C": {"PulseOffset": edit_pulse_1, "PulseDuration": 0.02},
        "D": {"PulseOffset": edit_pulse_off, "PulseDuration": 0.02}}

    meta_edited.set_dim_info(
        0,
        'DIM_EDIT',
        dim_info,
        dim_header)

    meta_edited.set_dim_info(1, 'DIM_DYN')
    meta_edited.set_standard_def("EditPulse", edit_pulse_val)

    return [gen_nifti_mrs_hdr_ext(data_short_te, dwelltime, meta_short_te, orientation, no_conj=True),
            gen_nifti_mrs_hdr_ext(data_edited, dwelltime, meta_edited, orientation, no_conj=True)],\
           [fout_str + '_hyper_short_te', fout_str + '_hyper_edited']
