"""spec2nii module containing functions specific to interpreting Bruker formats
Author: Tomas Psorn <tomaspsorn@isibrno.cz>
Copyright (C) 2021 Institute of Scientific Instruments of the CAS, v. v. i.
"""

from brukerapi.dataset import Dataset
from brukerapi.folders import Folder, TypeFilter
from brukerapi.splitters import SlicePackageSplitter, FrameGroupSplitter
from brukerapi.mergers import FrameGroupMerger
from brukerapi.exceptions import FilterEvalFalse
import numpy as np
import os
import sys
import pkg_resources
from spec2nii.nifti_orientation import NIFTIOrient
from spec2nii import nifti_mrs
import numpy as np

def read_bruker(args):
    """

    :param args:
    :return list imageOut: 
    :return list fileoutNames:
    """
    imageOut = []
    fileoutNames = []

    # for all Bruker datasets compliant all queries
    for data, properties in yield_bruker(args):
            orientation = NIFTIOrient(np.reshape(np.array(properties['affine']), (4,4)))
            imageOut.append(
                nifti_mrs.NIfTI_MRS(data,
                                   orientation.Q44,
                                   properties['dwell_s'],
                                   nifti_mrs.hdr_ext(
                                       properties['SpectrometerFrequency'],
                                       properties['ResonantNucleus']
                                   ))
            )           
            fileoutNames.append(properties['id'])

    return imageOut, fileoutNames

def yield_bruker(args):
    """

    If the path spectified by args.file is: 

    1/ Bruker dataset file (2dseq) - function yields its data and properties of the dataset
    2/ Directory - function yields data and properties and data of all datasets compliant to the queries 

    """
    # get location of the spec2nii Bruker properties configuration file
    bruker_properties_path = pkg_resources.resource_filename('spec2nii', 'bruker_properties.json')
    
    # get a list of queries to filter datasets
    queries = _get_queries(args)

    # case of Bruker dataset
    if os.path.isfile(args.file):
        d = Dataset(args.file, property_files=[bruker_properties_path])
        try:
            d.query(queries)
        except FilterEvalFalse:
            raise ValueError(f'Bruker dataset {d.path} is not suitable for conversion to mrs_nifti')
        yield from _proc_dataset(d)

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
                yield from _proc_dataset(d)

def _get_queries(args):
    """
    Returns a list of queries for filtering out only spectroscopic 2dseq datasets with a complex frame group

    """
    if args.mode == '2DSEQ':
        queries = ["@type=='2dseq'", "@is_spectroscopy==True", "@is_complex==True"]
    elif args.mode == 'FID':
        queries = ["@type=='fid'", "@is_spectroscopy==True"]
    return queries + args.query

def _proc_dataset(d):
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

    # some Bruker datasets do not have affine property
    if d.type == 'fid': if not 'affine' in properties: properties.update({'affine':np.identity(4)})
    
    yield data, properties

def _prep_data_svs(d):
    """
    Push the spectral dimension of the data array to the 3rd position for SVS data

    It is possible to use tuple as an axis argument of the expand_dims function since numpy>=1.18.0,
    we decided to use this triple call to avoid limiting numpy versions

    """
    data = np.expand_dims(d.data, axis=0)
    data = np.expand_dims(data, axis=0)
    data = np.expand_dims(data, axis=0)
    return data

def _prep_data_mrsi(d):
    """
    Push the spectral dimension of the data array to the 3rd position for CSI data

    """
    # push the spectral dimension to possition 2
    data = np.moveaxis(d.data, 0, 2)
    # add empty dimensions to push the spectral dimension to the 3rd index
    data = np.expand_dims(data, axis=2)
    return data