"""spec2nii module containing functions specific to interpreting Philips formats
Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford
"""
from datetime import datetime
from ast import literal_eval

import numpy as np

from spec2nii.nifti_orientation import NIFTIOrient, calc_affine
from spec2nii import nifti_mrs


def read_sdat_spar_pair(sdat_file, spar_file):

    spar_params = read_spar(spar_file)
    data = read_sdat(sdat_file,
                     spar_params['samples'],
                     spar_params['rows'])

    if data.ndim < 4:
        data = data.reshape((1, 1, 1) + data.shape)

    # Move to right handed frame
    data = data.conj()

    # Dwelltime
    dwelltime = 1.0 / float(spar_params["sample_frequency"])

    # Meta
    meta = spar_to_nmrs_hdrext(spar_params)
    meta.set_standard_def('OriginalFile', [sdat_file.name])

    # meta_used.set_dim_info(0, tag)

    # Orientation
    affine = _philips_orientation(spar_params)
    orientation = NIFTIOrient(affine)

    return [nifti_mrs.NIfTI_MRS(data, orientation.Q44, dwelltime, meta), ]


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
    raw_data = np.reshape(raw_data, (rows, samples)).squeeze()

    return raw_data


def _philips_orientation(params):
    '''Calculate the orientation affine from the spar parmaeters.'''

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
    obj = nifti_mrs.hdr_ext(cf,
                            spar_dict["nucleus"])

    # Some scanner information
    def set_standard_def(nifti_mrs_key, location, key, cast=None):
        try:
            if cast is not None:
                obj.set_standard_def(nifti_mrs_key, cast(getattr(location, key)))
            else:
                obj.set_standard_def(nifti_mrs_key, getattr(location, key))
        except AttributeError:
            pass

    obj.set_standard_def('Manufacturer', 'Philips')
    # set_standard_def('ManufacturersModelName',)
    # set_standard_def('DeviceSerialNumber', )
    set_standard_def('SoftwareVersions', spar_dict, 'equipment_sw_verions')

    # set_standard_def('InstitutionName', )
    # set_standard_def('InstitutionAddress', )

    # Some sequence information
    # obj.set_standard_def('SequenceName', )
    set_standard_def('ProtocolName', spar_dict, 'scan_id')

    # Some subject information
    try:
        obj.set_standard_def('PatientPosition', spar_dict['patient_position'] + ' ' + spar_dict['patient_orientation'])
    except AttributeError:
        pass
    set_standard_def('PatientName', spar_dict, 'patient_name')
    # set_standard_def('PatientWeight', , cast=float)
    set_standard_def('PatientDoB', spar_dict, 'patient_birth_date')
    # set_standard_def('PatientSex', )

    # Timing and sequence parameters
    obj.set_standard_def('EchoTime', float(spar_dict['echo_time']) * 1E-3)
    if spar_dict['spectrum_inversion_time'] > 0:
        obj.set_standard_def('InversionTime', spar_dict['spectrum_inversion_time'])
    # obj.set_standard_def('ExcitationFlipAngle', spar_dict['spectrum_inversion_time'])
    obj.set_standard_def('RepetitionTime', spar_dict['repetition_time'] / 1E3)
    # TO DO  - nibabel might need updating.
    obj.set_standard_def('TxOffset', spar_dict['offset_frequency'] / cf)

    # Conversion information
    obj.set_standard_def('ConversionMethod', 'spec2nii')
    conversion_time = datetime.now().isoformat(sep='T', timespec='milliseconds')
    obj.set_standard_def('ConversionTime', conversion_time)

    return obj


# From VESPA - BSD license. Ask for licence text from Brian.
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
        elif expon == 0 and sign == 0:
            f.append(0)
        else:
            f.append(0)
            # may want to raise an exception here ...

    return f
