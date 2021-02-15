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
    imageOut = []
    fileoutNames = []
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

    :param args:
    :return:
    """
    # get location of the spec2nii Bruker properties configuration file
    bruker_properties_path = pkg_resources.resource_filename('spec2nii', 'bruker_properties.json')
    
    if os.path.isfile(args.file):
        d = Dataset(args.file, property_files=[bruker_properties_path])
        yield from _proc_dataset(d)
    elif os.path.isdir(args.file):
        # select datasets to convert
        queries = _get_queries(args)

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
                print("doing {}".format(str(dataset.path)))
                yield from _proc_dataset(d)

def _get_queries(args):
    """
    returns a list of queries for filtering out only spectroscopic 2dseq datasets with a complex frame group

    """
    queries = ["@type=='2dseq'", "@is_spectroscopy==True", "@is_complex==True"]
    return queries + args.query

def _proc_dataset(d):

    if d.num_slice_packages > 1:
        for d_ in SlicePackageSplitter().split(d):
            yield from _proc_dataset(d_)
    else:
        if '<FG_COMPLEX>' in d.dim_type:
            d = FrameGroupMerger().merge(d, 'FG_COMPLEX')
        # prepare data based on type of experiment
        if 'spectroscopic' in d.dim_type[0]:
            if 'spatial' in d.dim_type:
                data = _prep_data_mrsi(d)
            else:
                data = _prep_data_svs(d)
        else:
            data = d.data

        yield data, d.to_dict()

def _prep_data_svs(d):
    """
    Add empty dimensions to push the temporal dimension to the 3rd position

    It is possible to use tuple as an axis argument of the expand_dims function since numpy>=1.18.0,
    we decided to use this triple call to avoid limiting numpy versions

    :param d:
    :return:
    """
    data = np.expand_dims(d.data, axis=0)
    data = np.expand_dims(data, axis=0)
    data = np.expand_dims(data, axis=0)
    return data

def _prep_data_mrsi(d):
    # push the temporal dimension to possition 2
    data = np.moveaxis(d.data, 0, 2)
    # add empty dimensions to push the temporal dimension to the 3rd index
    data = np.expand_dims(data, axis=2)
    return data