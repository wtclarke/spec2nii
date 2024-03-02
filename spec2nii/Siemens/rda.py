""" spec2nii module containing functions specific to Siemens rda format
Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2022 University of Oxford
"""
import re
from datetime import datetime

import numpy as np

from nifti_mrs.create_nmrs import gen_nifti_mrs_hdr_ext
from nifti_mrs.hdr_ext import Hdr_Ext

from spec2nii.dcm2niiOrientation.orientationFuncs import dcm_to_nifti_orientation
from spec2nii import __version__ as spec2nii_ver


class MRSINotHandledError(Exception):
    pass


def _locale_float(x):
    """Handle locale specific floating point representations in header

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
            try:
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
            except UnicodeDecodeError:
                print('Trying latin-1 encoding.')
                if hdr_st.search(line.decode('latin-1')):
                    pass
                    # print('header found')
                elif hdr_end.search(line.decode('latin-1')):
                    # print('header end')
                    break
                else:
                    match = hdr_val.search(line.decode('latin-1'))
                    if len(match.groups()) < 2:
                        hdr[match[1]] = None
                    else:
                        hdr[match[1]] = match[2]
        if verbose:
            print(hdr)

        data = np.fromfile(fp)

    data_cmplx = data[0::2] + 1j * data[1::2]

    # MRSI
    if (int(hdr['CSIMatrixSize[0]'])
            * int(hdr['CSIMatrixSize[1]'])
            * int(hdr['CSIMatrixSize[2]']))\
            > 1:
        data_cmplx = data_cmplx.reshape((
            int(hdr['CSIMatrixSize[2]']),
            int(hdr['CSIMatrixSize[1]']),
            int(hdr['CSIMatrixSize[0]']),
            int(hdr['VectorSize'])))
        data_cmplx = np.moveaxis(data_cmplx, (0, 1, 2), (2, 1, 0))
        data_shape = data_cmplx.shape[:3]

        imagePositionPatient = np.asarray([
            _locale_float(hdr['PositionVector[0]']),
            _locale_float(hdr['PositionVector[1]']),
            _locale_float(hdr['PositionVector[2]'])])

        half_shift = True
    # SVS
    else:
        data_cmplx = data_cmplx.reshape((1, 1, 1) + data_cmplx.shape)
        data_shape = (1, 1, 1)

        try:
            imagePositionPatient = np.asarray([
                _locale_float(hdr['VOIPositionSag']),
                _locale_float(hdr['VOIPositionCor']),
                _locale_float(hdr['VOIPositionTra'])])
        except KeyError:
            imagePositionPatient = np.asarray([
                _locale_float(hdr['PositionVector[0]']),
                _locale_float(hdr['PositionVector[1]']),
                _locale_float(hdr['PositionVector[2]'])])

        half_shift = False

    dwelltime = _locale_float(hdr['DwellTime']) / 1E6

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
                                                    data_shape,
                                                    verbose=verbose,
                                                    half_shift=half_shift)

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
    def set_standard_def(nifti_mrs_key, key, cast=None):
        if key in hdr\
                and hdr[key] is not None\
                and hdr[key]:
            if cast is not None:
                obj.set_standard_def(nifti_mrs_key, cast(hdr[key]))
            else:
                obj.set_standard_def(nifti_mrs_key, hdr[key])

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
    set_standard_def('SequenceName', 'SequenceName')
    # 'ProtocolName'
    set_standard_def('ProtocolName', 'ProtocolName')
    # # 5.4 Sequence information
    # 'PatientPosition'
    set_standard_def('PatientPosition', 'PatientPosition')
    # 'PatientName'
    set_standard_def('PatientName', 'PatientName')
    # 'PatientID'
    set_standard_def('PatientID', 'PatientID')
    # 'PatientWeight'
    set_standard_def('PatientWeight', 'PatientWeight', cast=_locale_float)
    # 'PatientDoB'
    set_standard_def('PatientDoB', 'PatientBirthDate')
    # 'PatientSex'
    set_standard_def('PatientSex', 'PatientSex')

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
