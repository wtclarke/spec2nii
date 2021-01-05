"""spec2nii module containing functions specific to interpreting jmrui formats
Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford
"""
import numpy as np
import re
from spec2nii.nifti_orientation import NIFTIOrient
from os.path import basename, splitext
from spec2nii import nifti_mrs
from datetime import datetime


def jmrui_txt(args):
    """Main function for handling jmrui text files."""

    data, header = readjMRUItxt(args.file)

    newshape = (1, 1, 1) + data.shape
    data = data.reshape(newshape)

    # meta-data
    dwelltime = header['SamplingInterval'] * 1E-3

    meta = jmrui_hdr_to_obj(header)
    meta.set_standard_def('OriginalFile', [basename(args.file), ])

    # Read optional affine file
    if args.affine:
        affine = np.loadtxt(args.affine)
    else:
        tmp = np.array([10000, 10000, 10000, 1])
        affine = np.diag(tmp)

    nifti_orientation = NIFTIOrient(affine)

    img_out = [nifti_mrs.NIfTI_MRS(data,
                                   nifti_orientation.Q44,
                                   dwelltime,
                                   meta), ]

    # File names
    if args.fileout:
        fname_out = [args.fileout, ]
    else:
        fname_out = [splitext(basename(args.file))[0], ]

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
    with open(filename, 'r') as txtfile:
        for line in txtfile:
            headerComp = headerRe.match(line)
            if headerComp:
                value = headerComp[2].strip()
                header.update({headerComp[1]: num(value)})

            signalIndices = signalRe.match(line)
            if signalIndices:
                recordData = True
                continue

            if recordData:
                curr_data = line.split()
                if len(curr_data) > 2:
                    curr_data = curr_data[:2]
                data.append(list(map(float, curr_data)))

    # Reshape data
    data = np.concatenate([np.array(i) for i in data])
    data = (data[0::2] + 1j * data[1::2]).astype(np.complex)

    return data, header


def num(s):
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return s


def jmrui_hdr_to_obj(header):
    """Translate jMRUI txt header to NIfTI MRS"""

    nucleus = id_nucleus(header)

    meta = nifti_mrs.hdr_ext(header['TransmitterFrequency'],
                             nucleus)

    meta.set_standard_def('ManufacturersModelName', header['Spectrometer'])
    meta.set_standard_def('PatientName', header['NameOfPatient'])

    meta.set_standard_def('ConversionMethod', 'spec2nii')
    conversion_time = datetime.now().isoformat(sep='T', timespec='milliseconds')
    meta.set_standard_def('ConversionTime', conversion_time)

    meta.set_user_def(key='AdditionalInformation',
                      value=header['AdditionalInformation'],
                      doc='jMRUI AdditionalInformation field.')

    return meta


nuclei = ['1H', '2H', '3H', '13C', '23NA', '31P', '129XE']
gyro_ratios = [42.577, 6.536, 45.415, 10.7084, 11.262, 17.235, -11.777]


def id_nucleus(header):
    """Identify the nucleus string.
    TypeOfNucleus is a header key but
    some files seem to not have it in."""

    if 'TypeOfNucleus' in header:
        return header['TypeOfNucleus']
    else:
        gyrom_mag_r = header['TransmitterFrequency'] / header['MagneticField']
        for nn, gr in zip(nuclei, gyro_ratios):
            if np.isclose(gr, gyrom_mag_r, atol=5E0):
                return nn

        print('Cannot identify nucleus setting 1H default.')
        return '1H'
