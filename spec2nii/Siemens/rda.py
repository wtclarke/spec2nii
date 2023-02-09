""" spec2nii module containing functions specific to Siemens rda format
Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2022 University of Oxford
"""
import re
from datetime import datetime
import warnings

import numpy as np

from nifti_mrs.create_nmrs import gen_nifti_mrs_hdr_ext
from nifti_mrs.hdr_ext import Hdr_Ext

from spec2nii.dcm2niiOrientation.orientationFuncs import dcm_to_nifti_orientation
from spec2nii import __version__ as spec2nii_ver


class MRSINotHandledError(Exception):
    pass


def _locale_float(x):
    """Handle locale specific flaoting point representations in header

    :param x: Header value as string with either . or , decimal separator
    :type x: str
    :return: Value cast to float
    :rtype: float
    """
    try:
        return float(x)
    except ValueError as e:
        if ',' in x:
            return float(x.replace(',', '.'))
        else:
            raise e


def convert_rda(rda_path, fname_out, verbose):
    hdr_st = re.compile(r'>>> Begin of header <<<')
    hdr_val = re.compile(r'^([\d\w\[\],]+): (.*)\r\n')
    hdr_end = re.compile(r'>>> End of header <<<')
    hdr = {}

    with open(rda_path, 'rb') as fp:
        for line in fp:
            if hdr_st.search(line.decode()):
                pass
                # print('header found')
            elif hdr_end.search(line.decode()):
                # print('header end')
                break
            else:
                match = hdr_val.search(line.decode())
                if len(match.groups()) < 2:
                    hdr[match[1]] = None
                else:
                    hdr[match[1]] = match[2]
        if verbose:
            print(hdr)

        data = np.fromfile(fp)

    if (int(hdr['CSIMatrixSize[0]'])
            * int(hdr['CSIMatrixSize[1]'])
            * int(hdr['CSIMatrixSize[2]']))\
            > 1:
        raise MRSINotHandledError('MRSI is currently not handled in the RDA format. Test data needed.')

    data_cmplx = data[0::2] + 1j * data[1::2]
    data_cmplx = data_cmplx.reshape((1, 1, 1) + data_cmplx.shape)
    dwelltime = _locale_float(hdr['DwellTime']) / 1E6

    warnings.warn(
        'The orientation calculations for rda data is mostly untested.'
        ' Please contribute test data if you can!')

    imagePositionPatient = np.asarray([
        _locale_float(hdr['PositionVector[0]']),
        _locale_float(hdr['PositionVector[1]']),
        _locale_float(hdr['PositionVector[2]'])])

    imageOrientationPatient = np.asarray([
        [_locale_float(hdr['RowVector[0]']), _locale_float(hdr['ColumnVector[0]'])],
        [_locale_float(hdr['RowVector[1]']), _locale_float(hdr['ColumnVector[1]'])],
        [_locale_float(hdr['RowVector[2]']), _locale_float(hdr['ColumnVector[2]'])]]).T

    xyzMM = np.asarray([
        _locale_float(hdr['PixelSpacingRow']),
        _locale_float(hdr['PixelSpacingCol']),
        _locale_float(hdr['SliceThickness'])])

    currNiftiOrientation = dcm_to_nifti_orientation(imageOrientationPatient,
                                                    imagePositionPatient,
                                                    xyzMM,
                                                    (1, 1, 1),
                                                    verbose=verbose)

    meta = extractRdaMetadata(hdr)
    meta.set_standard_def('OriginalFile', [rda_path.name])

    if fname_out:
        name = fname_out
    else:
        name = rda_path.stem

    return [gen_nifti_mrs_hdr_ext(data_cmplx, dwelltime, meta, currNiftiOrientation.Q44, no_conj=True), ], [name, ]


def extractRdaMetadata(hdr):
    """Extract information from the rda header dict to insert into the json header ext.

    :param hdr: Header read from rda file
    :type hdr: dict
    :return: Complete nifti_mrs header object
    :rtype: nifti_mrs.nifti_mrs.Hdr_Ext
    """

    # Extract required metadata and create hdr_ext object
    obj = Hdr_Ext(
        _locale_float(hdr['MRFrequency']),
        hdr['Nucleus'])

    # Standard defined metadata
    def set_standard_def(nifti_mrs_key, location, key, cast=None):
        if key in location\
                and location[key] is not None:
            if cast is not None:
                obj.set_standard_def(nifti_mrs_key, cast(location[key]))
            else:
                obj.set_standard_def(nifti_mrs_key, location[key])

    # # 5.1 MRS specific Tags
    # 'EchoTime'
    obj.set_standard_def('EchoTime', _locale_float(hdr['TE']) * 1E-3)
    # 'RepetitionTime'
    obj.set_standard_def('RepetitionTime', _locale_float(hdr['TR']) * 1E-3)
    # 'InversionTime'
    obj.set_standard_def('InversionTime', _locale_float(hdr['TI']) * 1E-3)
    # 'MixingTime'
    obj.set_standard_def('MixingTime', _locale_float(hdr['TM']) * 1E-3)
    # 'ExcitationFlipAngle'
    obj.set_standard_def('ExcitationFlipAngle', _locale_float(hdr['FlipAngle']))
    # 'TxOffset'
    # 'VOI'
    # 'WaterSuppressed'
    # 'WaterSuppressionType'
    # 'SequenceTriggered'

    # # 5.2 Scanner information
    # 'Manufacturer'
    obj.set_standard_def('Manufacturer', 'Siemens')
    # 'ManufacturersModelName'
    obj.set_standard_def('ManufacturersModelName', hdr['ModelName'])
    # 'DeviceSerialNumber'
    obj.set_standard_def('DeviceSerialNumber', hdr['DeviceSerialNumber'])
    # 'SoftwareVersions'
    obj.set_standard_def('SoftwareVersions', hdr['SoftwareVersion[0]'])
    # 'InstitutionName'
    obj.set_standard_def('InstitutionName', hdr['InstitutionName'])
    # 'InstitutionAddress'
    # 'TxCoil'
    obj.set_standard_def('TxCoil', hdr['TransmitCoil'])
    # 'RxCoil'

    # # 5.3 Sequence information
    # 'SequenceName'
    obj.set_standard_def('SequenceName', hdr['SequenceName'])
    # 'ProtocolName'
    obj.set_standard_def('ProtocolName', hdr['ProtocolName'])
    # # 5.4 Sequence information
    # 'PatientPosition'
    obj.set_standard_def('PatientPosition', hdr['PatientPosition'])
    # 'PatientName'
    obj.set_standard_def('PatientName', hdr['PatientName'])
    # 'PatientID'
    obj.set_standard_def('PatientID', hdr['PatientID'])
    # 'PatientWeight'
    obj.set_standard_def('PatientWeight', _locale_float(hdr['PatientWeight']))
    # 'PatientDoB'
    obj.set_standard_def('PatientDoB', hdr['PatientBirthDate'])
    # 'PatientSex'
    obj.set_standard_def('PatientSex', hdr['PatientSex'])

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
