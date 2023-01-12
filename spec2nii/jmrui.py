"""spec2nii module containing functions specific to interpreting jmrui formats
Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford
"""
import re
from datetime import datetime

import numpy as np

from nifti_mrs.create_nmrs import gen_nifti_mrs_hdr_ext
from nifti_mrs.hdr_ext import Hdr_Ext

from spec2nii.nifti_orientation import NIFTIOrient
from spec2nii import __version__ as spec2nii_ver


default_nuc_index = [None, "1H", "31P", "13C", "19F", "23NA"]
nuclei = ['1H', '2H', '3H', '13C', '23NA', '31P', '129XE']
gyro_ratios = [42.577, 6.536, 45.415, 10.7084, 11.262, 17.235, -11.777]


def id_nucleus(nucleus_hdr, trans_freq, b0):
    """Identify the nucleus string. TypeOfNucleus is a header key but
    some files seem to not have it in and sometimes it is an int

    :param: nucleus_hdr value of the TypeOfNucleus header key
    :param trans_freq: Transmitter/central frequency
    :param b0: Field strength in tesla.
    """

    if nucleus_hdr is not None\
            and isinstance(nucleus_hdr, str):
        return nucleus_hdr
    elif nucleus_hdr is not None\
            and isinstance(nucleus_hdr, (int, float))\
            and int(nucleus_hdr) > -1\
            and default_nuc_index[int(nucleus_hdr)] is not None:
        return default_nuc_index[int(nucleus_hdr)]
    else:
        try:
            gyrom_mag_r = trans_freq / b0
        except ZeroDivisionError:
            print('Cannot identify nucleus, using 1H default.')
            return '1H'

        if gyrom_mag_r > 1E4:
            gyrom_mag_r /= 1E6
        for nn, gr in zip(nuclei, gyro_ratios):
            if np.isclose(gr, gyrom_mag_r, atol=5E0):
                return nn

        print('Cannot identify nucleus, using 1H default.')
        return '1H'


def jmrui_format(args):
    '''Identify if text or binary format jmrui.'''
    if args.file.suffix == '.txt':
        return jmrui_txt(args)
    elif args.file.suffix == '.mrui':
        return jmrui_mrui(args)
    else:
        raise ValueError('jMRUI file format must be .txt or .mrui format.')


def jmrui_mrui(args):
    '''Process .mrui format files.'''

    data, header, str_info = read_mrui(args.file)

    newshape = (1, 1, 1) + data.shape
    data = data.reshape(newshape)

    # meta-data
    dwelltime = header['sampling_interval'] * 1E-3

    meta = jmrui_hdr_to_obj_mrui(header, str_info)
    meta.set_standard_def('OriginalFile', [args.file.name, ])

    if data.ndim > 4:
        meta.set_dim_info(0, 'DIM_USER_0', info='jMRUI frames')

    # Read optional affine file
    if args.affine:
        affine = np.loadtxt(args.affine)
    else:
        tmp = np.array([10000, 10000, 10000, 1])
        affine = np.diag(tmp)

    nifti_orientation = NIFTIOrient(affine)

    img_out = [gen_nifti_mrs_hdr_ext(
        data,
        dwelltime,
        meta,
        nifti_orientation.Q44,
        no_conj=True), ]

    # File names
    if args.fileout:
        fname_out = [args.fileout, ]
    else:
        fname_out = [args.file.stem, ]

    # Place in data output format
    return img_out, fname_out


def read_mrui(file_path):
    '''Read the header and data from a .mrui format file.

    :param Path file: file_path pathlib Path object to .mrui file.

    :return: Data in numpy format (npoints x nframes)
    :return: header dict
    :return: Additional information sorted as strings in file
    '''
    with open(file_path, 'br') as fp:
        # Header and data is stored as big endian double
        dt = np.dtype(np.float64)
        dt = dt.newbyteorder('>')

        # Read header values
        hdr = np.fromfile(fp, dtype=dt, count=13, offset=0)

        # Read data
        fp.seek(0, 2)
        file_size = fp.tell()
        n_frames = int(np.floor((file_size - 512) / (8 * hdr[1] * 2)))

        fp.seek(512, 0)
        data_count = int(n_frames * hdr[1] * 2)
        data = np.fromfile(fp, dtype=dt, count=data_count, offset=0)

        # Read final string info
        file_str = fp.read().decode("utf-8")

    # Sort header information
    header = {
        'type_of_sig': hdr[0],
        'number_of_points': int(hdr[1]),
        'sampling_interval': hdr[2],
        'begin_time ': hdr[3],
        'zero_order_phs': hdr[4],
        'transmitter_frequency': hdr[5],
        'magnetic_field': hdr[6],
        'type_of_nucleus': hdr[7],
        'reference_frequency_hz': hdr[8],
        'reference_frequency_ppm': hdr[9],
        'fid_or_echo': hdr[10],
        'apodizing': hdr[11],
        'num_zeros_view': hdr[12],
    }

    data = data[0::2] + 1j * data[1::2]

    data = data.reshape((n_frames, header['number_of_points'])).T
    data = data.squeeze()
    data = data.conj()

    return data, header, file_str


def jmrui_hdr_to_obj_mrui(header, str_info):
    """Translate jMRUI mrui header to NIfTI MRS"""

    nucleus = id_nucleus(header['type_of_nucleus'],
                         header['transmitter_frequency'],
                         header['magnetic_field'])

    meta = Hdr_Ext(
        header['transmitter_frequency'],
        nucleus)

    # meta.set_standard_def('ManufacturersModelName', header['Spectrometer'])
    # meta.set_standard_def('PatientName', header['NameOfPatient'])
    meta.set_standard_def('TxOffset', header['reference_frequency_ppm'])

    meta.set_standard_def('ConversionMethod', 'spec2nii')
    conversion_time = datetime.now().isoformat(sep='T', timespec='milliseconds')
    meta.set_standard_def('ConversionTime', conversion_time)

    meta.set_user_def(key='AdditionalInformation',
                      value=str_info,
                      doc='jMRUI .mrui string field.')

    return meta


def jmrui_txt(args):
    """Function for handling jmrui text files."""

    data, header = readjMRUItxt(args.file)

    newshape = (1, 1, 1) + data.shape
    data = data.reshape(newshape)

    # meta-data
    dwelltime = header['SamplingInterval'] * 1E-3

    meta = jmrui_hdr_to_obj(header)
    meta.set_standard_def('OriginalFile', [args.file.name, ])

    if data.ndim > 4:
        meta.set_dim_info(0, 'DIM_USER_0', info='jMRUI signals')

    # Read optional affine file
    if args.affine:
        affine = np.loadtxt(args.affine)
    else:
        tmp = np.array([10000, 10000, 10000, 1])
        affine = np.diag(tmp)

    nifti_orientation = NIFTIOrient(affine)

    img_out = [gen_nifti_mrs_hdr_ext(
        data,
        dwelltime,
        meta,
        nifti_orientation.Q44,
        no_conj=True), ]

    # File names
    if args.fileout:
        fname_out = [args.fileout, ]
    else:
        fname_out = [args.file.stem, ]

    # Place in data output format
    return img_out, fname_out


def readjMRUItxt(filename):
    """
    Read .txt format file
    Parameters
    ----------
    filename : string
        Name of jmrui .txt file

    Returns
    -------
    array-like
        Complex data
    dict
        Header information

    """
    signalRe = re.compile(r'Signal (\d{1,}) out of (\d{1,}) in file')
    headerRe = re.compile(r'(\w*):(.*)')
    header = {}
    data   = []
    recordData = False
    signal_index = 0
    with open(filename, 'r') as txtfile:
        for line in txtfile:
            headerComp = headerRe.match(line)
            if headerComp:
                value = headerComp[2].strip()
                header.update({headerComp[1]: num(value)})

            signalIndices = signalRe.match(line)
            if signalIndices:
                recordData = True
                signal_index += 1
                continue

            if recordData:
                curr_data = line.split()
                if len(curr_data) > 2:
                    curr_data = curr_data[:2]
                data.append(list(map(float, curr_data)))

    # Reshape data
    data = np.concatenate([np.array(i) for i in data])
    data = (data[0::2] + 1j * data[1::2]).astype(complex)
    data = data.reshape((signal_index, -1)).T.squeeze()
    data = data.conj()

    return data, header


def num(s):
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return s


def check_transmitter_frequency_in_header(header):
    """Check that the central frequency in the header is in MHz.
    Correct if not.

    :param header: Header dict containing 'TransmitterFrequency' key.
    :type header: dict
    :return: Checked and modified header
    :rtype: dict
    """
    if header['TransmitterFrequency'] > 1E5:
        header['TransmitterFrequency'] /= 1E6

    return header


def jmrui_hdr_to_obj(header):
    """Translate jMRUI txt header to NIfTI MRS"""

    header = check_transmitter_frequency_in_header(header)

    if 'TypeOfNucleus' in header:
        nucleus = id_nucleus(header['TypeOfNucleus'],
                             header['TransmitterFrequency'],
                             header['MagneticField'])
    else:
        nucleus = id_nucleus(None,
                             header['TransmitterFrequency'],
                             header['MagneticField'])

    meta = Hdr_Ext(
        float(header['TransmitterFrequency']),
        nucleus)

    if 'Spectrometer' in header:
        meta.set_standard_def('ManufacturersModelName', header['Spectrometer'])
    if 'NameOfPatient' in header:
        meta.set_standard_def('PatientName', header['NameOfPatient'])

    meta.set_standard_def('ConversionMethod', f'spec2nii v{spec2nii_ver}')
    conversion_time = datetime.now().isoformat(sep='T', timespec='milliseconds')
    meta.set_standard_def('ConversionTime', conversion_time)

    if 'AdditionalInformation' in header:
        meta.set_user_def(key='AdditionalInformation',
                          value=header['AdditionalInformation'],
                          doc='jMRUI AdditionalInformation field.')
    if 'SignalNames' in header:
        meta.set_user_def(key='SignalNames',
                          value=header['SignalNames'],
                          doc='jMRUI SignalNames field.')

    return meta
