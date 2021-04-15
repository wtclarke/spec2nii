"""spec2nii module wrapping an importer for the Varian file format

"Author": Jack J. Miller, jack.miller@physics.org
Copyright (C) 2021 University of Oxford and University of Aarhus
"""

import numpy as np
from spec2nii.nifti_orientation import NIFTIOrient
from spec2nii import nifti_mrs
from datetime import datetime
from os.path import basename, splitext
from spec2nii import __version__ as spec2nii_ver
from datetime import datetime

import re
import spec2nii.varian as v


def read_varian(args):
    """read_varian -- load a varian fid. Note that this format is very flexible and some rather large assumptions are made: at the moment, this assumes little.
    :param file: path to the varian .fid directory, containing fid and procpar files.
    returns img_out, fnemae_out
    """
    dic, data = v.read(args.file)

    # extract number of coils
    number_of_coils = 0
    for i in dic['procpar']['rcvrs']['values'][0]:
        if re.search(i,'y'):
            number_of_coils += 1

    # number of time points -- probably
    number_of_time_points = float(dic['procpar']['arraydim']['values'][0]) / number_of_coils
    if (not number_of_time_points.is_integer()):
        raise ValueError('Coil reshaping failed')

    number_of_time_points = int(number_of_time_points)

    # spectral number of points
    number_of_spectral_points = int(int(dic['np'])/2);

    # reshape
    newshape = (1,1,1, number_of_spectral_points, number_of_coils, number_of_time_points)
    data = data.transpose();
    data = np.resize(data, newshape);

    # extract additional spectral metadata
    dwelltime = 1.0 / float(dic['procpar']['sw']['values'][0])
    imagingfreq=float(dic['procpar']['sfrq']['values'][0])
    nucleus = dic['procpar']['tn']['values'][0]
    nucleus = nucleus[1:]+ nucleus[0] #reshape to be in the form '13C' rather than 'C13'


    try:
        repitition_time = float(dic['procpar']['tr']['values'][0])
    except KeyError:
        pass

    try:
        echotime = float(dic['procpar']['te']['values'][0])#In ms if 'tis there
    except KeyError:
        echotime = float(dic['procpar']['pw']['values'][0]) + float(dic['procpar']['alfa']['values'][0])

    # extract voxel positions if available
    # NB:  0d spect -- thk etc not defined
    #      1d slice -- thk, and then pss along the usual (theta, phi, psi) angles (and orient)
    #      SVS voxel - vorient, vpsi, vphi, vtheta: define angle of voxel;
    #                - pos1, pos2, pos3           : define position of voxel
    #                - vox1, vox2, vox3, thkunit  : define size of voxel (no thkunit is mm)

    #TODO: Jack should implement the affine matrix correctly -- see philips example

    affine = np.diag(np.array([10000, 10000, 10000, 1])) # 10 m slab for now....
    orientation = NIFTIOrient(affine)

    # create object
    meta = nifti_mrs.hdr_ext(imagingfreq, nucleus)
    meta.set_standard_def('ConversionMethod', f'spec2nii v{spec2nii_ver}')
    meta.set_standard_def('EchoTime', echotime);
    meta.set_standard_def('RepetitionTime', repitition_time)
    meta.set_standard_def('Manufacturer', 'Varian')
    meta.set_standard_def('ProtocolName', dic['procpar']['seqfil']['values'][0])
    meta.set_standard_def('PatientName', dic['procpar']['comment']['values'][0])

    conversion_time = datetime.now().isoformat(sep='T', timespec='milliseconds')
    meta.set_standard_def('ConversionTime', conversion_time)
    meta.set_standard_def('OriginalFile', [basename(args.file)])

    # k-space
    meta.set_standard_def('kSpace', [False, False, False])

    # set tag dimensions
    meta.set_dim_info(0, "DIM_COIL")
    meta.set_dim_info(1, "DIM_DYN")

    # File names
    if args.fileout:
        fname_out = [args.fileout, ]
    else:
        fname_out = [splitext(basename(args.file))[0], ]

    return [ nifti_mrs.NIfTI_MRS(data, orientation.Q44, dwelltime, meta), ], fname_out
