"""spec2nii module containing functions specific to interpreting Bruker formats
Dependent on the brukerapi package developed by Tomas Psorn.
https://github.com/isi-nmr/brukerapi-python

Author: Tomas Psorn <tomaspsorn@isibrno.cz>
        Will Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2021 Institute of Scientific Instruments of the CAS, v. v. i.
"""
import os
import pkg_resources
import warnings
from datetime import datetime

import numpy as np

from brukerapi.dataset import Dataset
from brukerapi.folders import Folder
from brukerapi.mergers import FrameGroupMerger
from brukerapi.exceptions import FilterEvalFalse

from nifti_mrs.create_nmrs import gen_nifti_mrs_hdr_ext
from nifti_mrs.hdr_ext import Hdr_Ext

from spec2nii.nifti_orientation import NIFTIOrient
from spec2nii import __version__ as spec2nii_ver

# Default dimension assignments.
fid_dimension_defaults = {
    'repetition': "DIM_DYN",
    'channel': "DIM_COIL"}


def read_bruker(args):
    """

    :param args:
    :return list imageOut:
    :return list fileoutNames:
    """
    imageOut = []
    fileoutNames = []

    # for all Bruker datasets compliant all queries
    for data, orientation, dwelltime, meta, name in yield_bruker(args):
        imageOut.append(
            gen_nifti_mrs_hdr_ext(
                data,
                dwelltime,
                meta,
                orientation.Q44,
                no_conj=True)
        )
        fileoutNames.append(name)

    return imageOut, fileoutNames


def yield_bruker(args):
    """

    If the path specified by args.file is:

    1/ Bruker dataset file (2dseq) - function yields its data and properties of the dataset
    2/ Directory - function yields data and properties and data of all datasets compliant to the queries

    """
    # get location of the spec2nii Bruker properties configuration file
    bruker_properties_path = pkg_resources.resource_filename('spec2nii', 'bruker_properties.json')
    bruker_fid_override_path = pkg_resources.resource_filename('spec2nii', 'bruker_fid_override.json')

    # get a list of queries to filter datasets
    queries = _get_queries(args)

    # case of Bruker dataset
    if os.path.isfile(args.file):
        d = Dataset(
            args.file,
            property_files=[bruker_fid_override_path, bruker_properties_path],
            parameter_files=['method'])
        try:
            d.query(queries)
        except FilterEvalFalse:
            raise ValueError(f'Bruker dataset {d.path} is not suitable for conversion to mrs_nifti')
        yield from _proc_dataset(d, args)

    # case of folder containing Bruker datasets
    elif os.path.isdir(args.file):

        # process individual datasets
        for dataset in Folder(args.file, dataset_state={
            "parameter_files": ['method'],
            "property_files": [bruker_properties_path]
        }).get_dataset_list_rec():
            with dataset as d:
                try:
                    d.query(queries)
                except FilterEvalFalse:
                    continue
                yield from _proc_dataset(d, args)


def _get_queries(args):
    """
    Returns a list of queries for filtering out only spectroscopic 2dseq datasets with a complex frame group

    """
    if args.mode == '2DSEQ':
        queries = ["@type=='2dseq'", "@is_spectroscopy==True", "@is_complex==True"]
    elif args.mode == 'FID':
        queries = ["@type=='fid'", "@is_spectroscopy==True"]
    return queries + args.query


def _proc_dataset(d, args):
    """
    Yield data and properties of a single dataset

    """
    # merge 2dseq complex frame group if present
    if d.is_complex and d.type == '2dseq':
        d = FrameGroupMerger().merge(d, 'FG_COMPLEX')

    # prepare the data array
    if d.is_svs:
        data = _prep_data_svs(d)
    elif d.is_mrsi:
        data = _prep_data_mrsi(d)
    else:
        data = d.data

    # get properties
    properties = d.to_dict()

    # Orientation information
    if d.type == 'fid':
        orientation = NIFTIOrient(_fid_affine_from_params(d))
    else:
        orientation = NIFTIOrient(np.reshape(np.array(properties['affine']), (4, 4)))

    # Meta data
    if d.type == 'fid':
        meta = _fid_meta(d, dump=args.dump_headers)
    else:
        meta = _2dseq_meta(d, dump=args.dump_headers)

    # Dwelltime - to do resolve this factor of 2 issue
    if d.type == 'fid':
        dwelltime = d.dwell_s * 2
    else:
        dwelltime = d.dwell_s * 2

    if args.fileout:
        name = args.fileout + '_' + d.id.rstrip('_')
    else:
        name = d.id.rstrip('_')

    yield data, orientation, dwelltime, meta, name


def _prep_data_svs(d):
    """
    Push the spectral dimension of the data array to the 3rd position for SVS data

    It is possible to use tuple as an axis argument of the expand_dims function since numpy>=1.18.0,
    we decided to use this triple call to avoid limiting numpy versions

    """
    data = d.data
    if d.type == 'fid':
        # Remove points acquired before echo
        data = data[d.points_prior_to_echo:, ...]

        # fid data appears to need to be conjugated for NIFTI-MRS convention
        data = data.conj()

    data = np.expand_dims(data, axis=0)
    data = np.expand_dims(data, axis=0)
    data = np.expand_dims(data, axis=0)
    return data


def _prep_data_mrsi(d):
    """
    Push the spectral dimension of the data array to the 3rd position for CSI data

    """
    data = d.data
    if d.type == 'fid':
        # Remove points acquired before echo
        data = data[d.points_prior_to_echo:, ...]

        # fid data appears to need to be conjugated for NIFTI-MRS convention
        data = data.conj()

    # push the spectral dimension to position 2
    data = np.moveaxis(data, 0, 2)
    # add empty dimensions to push the spectral dimension to the 3rd index
    data = np.expand_dims(data, axis=2)
    return data


def _fid_affine_from_params(d):
    """ First attempt to create 4x4 affine from fid headers"""
    warnings.warn('The orientation of bruker fid data is mostly untested.')

    orientation = np.squeeze(d.parameters['method']['PVM_VoxArrGradOrient'].value)
    shift = np.squeeze(d.parameters['method']['PVM_VoxArrPosition'].value)
    # csshift = np.squeeze(d.parameters['method']['PVM_VoxArrCSDisplacement'].value)
    # shift += csshift
    size = np.squeeze(d.parameters['method']['PVM_VoxArrSize'].value)
    affine = np.zeros((4, 4))
    affine[3, 3] = 1
    reorder = [1, 2, 0]  # [0, 1, 2] *[0, 2, 1]* [1, 2, 0]
    affine[:3, :3] = orientation[reorder, :].T * size[reorder]
    affine[:3, 3] = shift[reorder]

    return affine


def _2dseq_meta(d, dump=False):
    """ Extract information from method and acqp file into hdr_ext.

    :param d: Dataset
    :return: NIfTI MRS hdr ext object.
    """

    # Extract required metadata and create hdr_ext object
    cf = d.SpectrometerFrequency
    obj = Hdr_Ext(
        cf,
        d.ResonantNucleus)

    # # 5.1 MRS specific Tags
    # 'EchoTime'
    if hasattr(d, 'TE'):
        obj.set_standard_def('EchoTime', float(d.TE * 1E-3))
    elif hasattr(d, 'method_TE'):
        obj.set_standard_def('EchoTime', float(d.method_TE * 1E-3))
    # 'RepetitionTime'
    if hasattr(d, 'TR'):
        obj.set_standard_def('RepetitionTime', float(d.TR / 1E3))
    elif hasattr(d, 'method_TR'):
        obj.set_standard_def('RepetitionTime', float(d.method_TR / 1E3))
    # 'InversionTime'
    # 'MixingTime'
    # 'ExcitationFlipAngle'
    # 'TxOffset'
    # Bit of a guess, not sure of units.
    obj.set_standard_def('TxOffset', float(d.working_offset[0]))
    # 'VOI'
    # 'WaterSuppressed'
    # No apparent parameter stored in the SPAR info.
    # 'WaterSuppressionType'
    # 'SequenceTriggered'
    # # 5.2 Scanner information
    # 'Manufacturer'
    obj.set_standard_def('Manufacturer', 'Bruker')
    # 'ManufacturersModelName'
    # 'DeviceSerialNumber'
    # 'SoftwareVersions'
    obj.set_standard_def('SoftwareVersions', d.PV_version)
    # 'InstitutionName'
    # 'InstitutionAddress'
    # 'TxCoil'
    # 'RxCoil'
    # # 5.3 Sequence information
    # 'SequenceName'
    obj.set_standard_def('SequenceName', d.method_desc)
    # 'ProtocolName'
    # # 5.4 Sequence information
    # 'PatientPosition'
    # 'PatientName'
    obj.set_standard_def('PatientName', d.subj_id)
    # 'PatientID'
    # 'PatientWeight'
    # 'PatientDoB'
    # 'PatientSex'
    # # 5.5 Provenance and conversion metadata
    # 'ConversionMethod'
    obj.set_standard_def('ConversionMethod', f'spec2nii v{spec2nii_ver}')
    # 'ConversionTime'
    conversion_time = datetime.now().isoformat(sep='T', timespec='milliseconds')
    obj.set_standard_def('ConversionTime', conversion_time)
    # 'OriginalFile'
    obj.set_standard_def('OriginalFile', [str(d.path), ])
    # # 5.6 Spatial information
    # 'kSpace'
    obj.set_standard_def('kSpace', [False, False, False])

    # Stuff full headers into user fields
    if dump:
        for hdr_file in d.parameters:
            obj.set_user_def(key=hdr_file,
                             doc=f'Bruker {hdr_file} file.',
                             value=d.parameters[hdr_file].to_dict())

    # Tags
    unknown_count = 0
    for ddx, dim in enumerate(d.dim_type[1:]):
        if dim in fid_dimension_defaults:
            obj.set_dim_info(ddx, fid_dimension_defaults[dim])
        else:
            obj.set_dim_info(ddx, f'DIM_USER_{unknown_count}')
            unknown_count += 1

    return obj


def _fid_meta(d, dump=False):
    """ Extract information from method and acqp file into hdr_ext.

    :param d: Dataset
    :return: NIfTI MRS hdr ext object.
    """

    # Extract required metadata and create hdr_ext object
    cf = d.SpectrometerFrequency
    obj = Hdr_Ext(
        cf,
        d.ResonantNucleus)

    # # 5.1 MRS specific Tags
    # 'EchoTime'
    obj.set_standard_def('EchoTime', float(d.TE * 1E-3))
    # 'RepetitionTime'
    obj.set_standard_def('RepetitionTime', float(d.TR / 1E3))
    # 'InversionTime'
    # 'MixingTime'
    # 'ExcitationFlipAngle'
    # 'TxOffset'
    # Bit of a guess, not sure of units.
    obj.set_standard_def('TxOffset', float(d.working_offset[0]))
    # 'VOI'
    # 'WaterSuppressed'
    # No apparent parameter stored in the SPAR info.
    # 'WaterSuppressionType'
    # 'SequenceTriggered'
    # # 5.2 Scanner information
    # 'Manufacturer'
    obj.set_standard_def('Manufacturer', 'Bruker')
    # 'ManufacturersModelName'
    # 'DeviceSerialNumber'
    # 'SoftwareVersions'
    obj.set_standard_def('SoftwareVersions', d.PV_version)
    # 'InstitutionName'
    # 'InstitutionAddress'
    # 'TxCoil'
    # 'RxCoil'
    # # 5.3 Sequence information
    # 'SequenceName'
    obj.set_standard_def('SequenceName', d.method_desc)
    # 'ProtocolName'
    # # 5.4 Sequence information
    # 'PatientPosition'
    # 'PatientName'
    obj.set_standard_def('PatientName', d.subj_id)
    # 'PatientID'
    # 'PatientWeight'
    # 'PatientDoB'
    # 'PatientSex'
    # # 5.5 Provenance and conversion metadata
    # 'ConversionMethod'
    obj.set_standard_def('ConversionMethod', f'spec2nii v{spec2nii_ver}')
    # 'ConversionTime'
    conversion_time = datetime.now().isoformat(sep='T', timespec='milliseconds')
    obj.set_standard_def('ConversionTime', conversion_time)
    # 'OriginalFile'
    obj.set_standard_def('OriginalFile', [str(d.path), ])
    # # 5.6 Spatial information
    # 'kSpace'
    obj.set_standard_def('kSpace', [False, False, False])

    # Stuff full headers into user fields
    if dump:
        for hdr_file in d.parameters:
            obj.set_user_def(key=hdr_file,
                             doc=f'Bruker {hdr_file} file.',
                             value=d.parameters[hdr_file].to_dict())

    # Tags
    unknown_count = 0
    for ddx, dim in enumerate(d.dim_type[1:]):
        if dim in fid_dimension_defaults:
            obj.set_dim_info(ddx, fid_dimension_defaults[dim])
        else:
            obj.set_dim_info(ddx, f'DIM_USER_{unknown_count}')
            unknown_count += 1

    return obj
