"""spec2nii module wrapping an importer for the Varian file format

"Author": Jack J. Miller, jack.miller@physics.org
Copyright (C) 2021 University of Oxford and University of Aarhus
"""

import re
import warnings
from datetime import datetime
from os.path import basename, splitext

import numpy as np

from nifti_mrs.create_nmrs import gen_nifti_mrs_hdr_ext
from nifti_mrs.hdr_ext import Hdr_Ext

from spec2nii import __version__ as spec2nii_ver
from spec2nii.nifti_orientation import NIFTIOrient, calc_affine
import spec2nii.varian as v


def read_varian(args):
    """read_varian -- load a varian fid.
    Note that this format is very flexible and some rather large assumptions are made.
    At the moment, this assumes little.
    :param file: path to the varian .fid directory, containing fid and procpar files.
    returns img_out, fname_out
    """
    dic, data = v.read(args.file)

    # extract number of coils
    number_of_coils = 0
    for i in dic['procpar']['rcvrs']['values'][0]:
        if re.search(i, 'y'):
            number_of_coils += 1

    # number of time points -- probably
    number_of_time_points = float(dic['procpar']['arraydim']['values'][0]) / number_of_coils
    if (not number_of_time_points.is_integer()):
        raise ValueError('Coil reshaping failed')

    number_of_time_points = int(number_of_time_points)

    # spectral number of points
    number_of_spectral_points = int(int(dic['np']) / 2)

    # number of slices
    sequence_name = dic['procpar']['seqfil']['values'][0]
    if (sequence_name.count('SSel') or sequence_name.count('sliceshim')):  # 1D localised slice
        pass  # to be created

    # reshape
    newshape = (1, 1, 1, number_of_spectral_points, number_of_coils, number_of_time_points)
    data = data.transpose()
    data = np.resize(data, newshape)

    # extract additional spectral metadata
    dwelltime = 1.0 / float(dic['procpar']['sw']['values'][0])
    imagingfreq = float(dic['procpar']['sfrq']['values'][0])
    nucleus = dic['procpar']['tn']['values'][0]
    # reshape to be in the form '13C' rather than 'C13'
    nucleus = nucleus[1:] + nucleus[0]

    try:
        repetition_time = float(dic['procpar']['tr']['values'][0])
    except KeyError:
        pass

    try:
        echotime = float(dic['procpar']['te']['values'][0])  # In ms if 'tis there
    except KeyError:
        pass
    try:
        echotime = float(dic['procpar']['pw']['values'][0]) + float(dic['procpar']['alfa']['values'][0])
    except KeyError:
        pass
    try:
        echotime = float(dic['procpar']['p1']['values'][0]) + float(dic['procpar']['alfa']['values'][0])
    except KeyError:
        pass

    # Parse 3D localisation
    if (sequence_name.count('press') or sequence_name.count('steam')):
        affine = _varian_orientation_3d(dic)
    else:
        affine = np.diag(np.array([10000, 10000, 10000, 1]))  # 10 m slab for now....

    # TODO: Jack should implement the affine matrix correctly for all sequences

    orientation = NIFTIOrient(affine)

    # create object
    meta = Hdr_Ext(imagingfreq, nucleus)
    meta.set_standard_def('ConversionMethod', f'spec2nii v{spec2nii_ver}')
    meta.set_standard_def('EchoTime', echotime)
    meta.set_standard_def('RepetitionTime', repetition_time)
    meta.set_standard_def('Manufacturer', 'Varian')
    meta.set_standard_def('ProtocolName', dic['procpar']['seqfil']['values'][0])
    meta.set_standard_def('PatientName', dic['procpar']['comment']['values'][0])

    # stuff that is nice to have:
    try:
        meta.set_standard_def('SoftwareVersions', dic['procpar']['parver']['values'][0])
        meta.set_standard_def('TxCoil', dic['procpar']['rfcoil']['values'][0])
        meta.set_standard_def('RxCoil', dic['procpar']['rfcoil']['values'][0])
    except KeyError:
        warnings.warn('Expected standard metadata keying failed')
    try:
        meta.set_standard_def('InversionTime', float(dic['procpar']['ti']['values'][0]))
    except KeyError:
        pass
    try:
        meta.set_standard_def('ExcitationFlipAngle', float(dic['procpar']['flip1']['values'][0]))
    except KeyError:
        pass
    conversion_time = datetime.now().isoformat(sep='T', timespec='milliseconds')
    meta.set_standard_def('ConversionTime', conversion_time)
    meta.set_standard_def('OriginalFile', [basename(args.file)])

    # k-space
    meta.set_standard_def('kSpace', [False, False, False])

    # set tag dimensions
    meta.set_dim_info(0, "DIM_COIL")
    meta.set_dim_info(1, args.tag6)

    # Stuff full headers into user fields
    if args.dump_headers:
        meta.set_user_def(key='VarianProcpar',
                          doc='Varian procpar metadata.',
                          value=dic['procpar'])

    # File names
    if args.fileout:
        fname_out = [args.fileout, ]
    else:
        fname_out = [splitext(basename(args.file))[0], ]

    return [gen_nifti_mrs_hdr_ext(data, dwelltime, meta, orientation.Q44, no_conj=True), ], fname_out


def _varian_orientation_1d(params):
    '''Calculate 1d slice orientation from parameters struct
     extract voxel positions if available
     NB:  0d spect -- thk etc not defined
          1d slice -- thk, and then pss along the usual (theta, phi, psi) angles (and orient)
          SVS voxel - vorient, vpsi, vphi, vtheta: define angle of voxel;
                    - pos1, pos2, pos3           : define position of voxel
                    - vox1, vox2, vox3, thkunit  : define size of voxel (no thkunit is mm)
    '''
    warnings.warn('Not yet implemented')


def _varian_orientation_3d(params):
    '''Calculate single voxel spectroscopy orientation from parameters struct'''

    angle_lr = params['procpar']['vpsi']['values'][0]
    angle_ap = params['procpar']['vtheta']['values'][0]
    angle_hf = params['procpar']['vphi']['values'][0]
    angles = [angle_lr, angle_ap, -angle_hf]

    dim_lr = params['procpar']['vox1']['values'][0]
    dim_ap = params['procpar']['vox2']['values'][0]
    dim_hf = params['procpar']['vox3']['values'][0]
    dimensions = [dim_lr, dim_ap, dim_hf] * 1e-2

    shift_lr = params['procpar']['pos1']['values'][0]
    shift_ap = params['procpar']['pos2']['values'][0]
    shift_hf = params['procpar']['pos3']['values'][0]
    shift = [-shift_lr, -shift_ap, shift_hf] * 1e-2

    return calc_affine(angles, dimensions, shift)
