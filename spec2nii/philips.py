"""spec2nii module containing functions specific to interpreting Philips formats
Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford
"""

from ast import literal_eval
import numpy as np
from spec2nii.nifti_orientation import NIFTIOrient, calc_affine
import warnings


def read_sdat_spar_pair(sdat_file, spar_file):

    spar_params = read_spar(spar_file)
    data = read_sdat(sdat_file,
                     spar_params['samples'],
                     spar_params['rows'])

    # Dwelltime
    dwelltime = 1.0 / float(spar_params["sample_frequency"])

    # Meta - which will be written to the JSON side car requires
    # a few parameters. Set these first. Then set the ones from
    # the spar file
    cf = float(spar_params["synthesizer_frequency"]) / 1E6
    meta = {'ImagingFrequency': cf, 'Dwelltime': dwelltime}
    meta.update(spar_params)

    # Orientation
    affine = _philips_orientation(spar_params)
    orientation = NIFTIOrient(affine)

    return data, orientation, dwelltime, meta


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
            except ValueError or SyntaxError:
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
    warnings.warn('Philips orientation not yet tested.')

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
