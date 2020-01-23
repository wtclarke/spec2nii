import nibabel as nib
import os.path as op
import numpy as np

def writeNii(filename,outdir,data,positionInfo,dwellTime):

    # Form full path
    fullfilepath = op.join(outdir,filename+'.nii.gz')

    # reshape to have three singleton dimensions - this will have to be modified to accept CSI data
    newshape = (1,1,1)+data.shape
    data = data.reshape(newshape)

    # Create new nifti image
    newobj = nib.nifti2.Nifti2Image(data,positionInfo.Q44)

    # Write new header
    #newobj = updateHeader(newobj,positionInfo,dwellTime)
    pixDim = newobj.header['pixdim']
    pixDim[4] = dwellTime
    newobj.header['pixdim'] = pixDim
    #import pdb; pdb.set_trace()
    # From nii obj and write    
    nib.save(newobj,fullfilepath)

    return None