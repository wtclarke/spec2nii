""" spec2nii module for writing output sidecar JSON files
Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford 
"""
import json
import os.path as op
import numpy as np

def writeJSON(filename,outdir,meta_info):
    """ Write spectroscopy meta-data as json
    
    Args:
        filename (str): Output file name, less extension
        outdir (str): Output directory        
        meta_info (dict): Meta info to write to json.
    """
    # Form full path
    fullfilepath = op.join(outdir,filename+'.json')

    # Clean up dict
    for key in meta_info:
        if isinstance(meta_info[key],np.ndarray):
            meta_info[key] = meta_info[key].tolist()

        # print(f'{key} is type: {type(meta_info[key])}.')

    with open(fullfilepath, 'w') as fp:
        json.dump(meta_info, fp)
    

