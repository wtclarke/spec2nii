""" spec2nii module for writing output NIFTI images
Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford 
"""

import nibabel as nib
import os.path as op
import numpy as np

def writeNii(filename,outdir,data,positionInfo,dwelltime):
    """ Write spectroscopy data as NIFTI
    
    Args:
        filename (str): Output file name less extension
        outdir (str): Output directory
        data (np.ndarray): Time domain data (excluding spatial dimensions)
        positionInfo (NIFTIOrient obj): Orientation information.
        dwelltime (float): Dwelltime in seconds

    """
    # Form full path
    fullfilepath = op.join(outdir,filename+'.nii.gz')

    # reshape to have three singleton dimensions - this will have to be modified to accept CSI data
    newshape = (1,1,1)+data.shape
    data = data.reshape(newshape)

    # Create new nifti image
    newobj = nib.nifti2.Nifti2Image(data,positionInfo.Q44)

    # Write new header
    pixDim = newobj.header['pixdim']
    pixDim[4] = dwelltime
    newobj.header['pixdim'] = pixDim
    
    # From nii obj and write    
    nib.save(newobj,fullfilepath)

    return None