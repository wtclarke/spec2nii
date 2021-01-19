'''spec2nii, functions to interpret GE P files.
Very hacky interpretation of GE p files focused on use for SVS MRS
Cobbled together from three packages: Spant, Gannet and pfile-tools
I have no documentation for or experience with GE data, so most of this is trial & error!

Author: William Clarke <william.clarke@ndcn.ox.ac.uk>

Copyright (C) 2020 University of Oxford
'''

import warnings
from datetime import datetime
from copy import deepcopy

import numpy as np

from spec2nii.nifti_orientation import NIFTIOrient, calc_affine
from spec2nii import nifti_mrs


class IsMRSI(Exception):
    pass


def read_p_file(filename):
    '''Read and interpret a GE p-file.
    :param filename: Path to file
    '''

    hdr = _read_hdr(filename)

    if (hdr['xcsi'] * hdr['ycsi'] * hdr['zcsi']) > 1:
        raise IsMRSI(f"Only svs is handled currently. x:{hdr['xcsi']},y:{hdr['ycsi']},z:{hdr['zcsi']}.")

    data, ref_data = _read_data(filename, hdr)

    data = data.reshape((1, 1, 1) + data.shape)
    if ref_data is not None:
        ref_data = ref_data.reshape((1, 1, 1) + ref_data.shape)

    affine = _calc_affine(hdr)
    orientation = NIFTIOrient(affine)

    dwelltime = 1 / hdr['spec_width']

    # Metadata
    warnings.warn('Assuming 1H data!')
    meta = nifti_mrs.hdr_ext(float(hdr["ps_mps_freq"]) / 1E7,
                             '1H')
    meta.set_standard_def('EchoTime', hdr['te'])
    meta.set_standard_def('RepetitionTime', hdr['tr'])

    meta.set_standard_def('ConversionMethod', 'spec2nii')
    conversion_time = datetime.now().isoformat(sep='T', timespec='milliseconds')
    meta.set_standard_def('ConversionTime', conversion_time)

    meta_ref = deepcopy(meta)

    meta.set_dim_info(0, 'DIM_COIL')
    meta.set_dim_info(1, 'DIM_DYN')

    meta_ref.set_dim_info(0, 'DIM_COIL')
    meta_ref.set_dim_info(1, 'DIM_DYN')

    return nifti_mrs.NIfTI_MRS(data, orientation.Q44, dwelltime, meta),\
        nifti_mrs.NIfTI_MRS(ref_data, orientation.Q44, dwelltime, meta_ref)


def _calc_affine(hdrs):
    '''Calculate the affine for a GE p-file.'''
    warnings.warn("Orientation information for GE files is incomplete. Orientation not handled!")

    angles = [0, 0, 0]  # TO DO!

    dimensions = hdrs['vox_size']

    shift = hdrs['vox_offset']

    return calc_affine(angles, dimensions, shift)


def _read_data(filename, hdrs):
    '''Read the data component of a p-file.'''
    frames = int(hdrs['nframes'])
    echoes = int(hdrs['nechoes'])
    fid_pnts = int(hdrs['frame_size'])
    ref_frames = int(hdrs['rhuser19'])

    with open(filename, 'rb') as f:
        f.seek(0, 2)
        file_size = f.tell()
        f.seek(0)
        pnts = int((file_size - hdrs['off_data']) / 4)
        print(pnts)
        f.seek(hdrs['off_data'])
        data = np.fromfile(f, dtype=np.int32, count=pnts)

    coils = 0
    for n in np.arange(0, 8, 2):
        if hdrs['rcv'][n] != 0 or hdrs['rcv'][n + 1] != 0:
            coils = coils + 1 + hdrs['rcv'][n + 1] - hdrs['rcv'][n]
    if coils == 0:
        coils = 1

    expt_pts = coils * (frames * echoes + echoes) * fid_pnts * 2

    if expt_pts != pnts:
        warnings.warn("Unexpected number of data points.")
        print(f"Expecting {pnts} points based on file size and {expt_pts} based on header.")
        print(f"Coils {coils}")
        print(f"nframes {frames}")
        print(f"nechoes {echoes}")
        print(f"frame_size {fid_pnts}")
        print(f"w_frames {ref_frames}")

    data = data.reshape((-1, 2))
    data = data[:, 0] + 1j * data[:, 1]

    dynamics = frames * echoes + echoes
    data = data.reshape((coils, dynamics, fid_pnts))

    # data = np.conj(data)

    # Remove empty frame
    data = data[:, 1:, ]

    # Make order fid_pnts, coils, dynamics
    data = np.moveaxis(data, (0, 1, 2), (1, 2, 0))
    if ref_frames == 0:
        ref_data = None
    else:
        ref_data = data[:, :, 0:ref_frames].copy()
    data = data[:, :, ref_frames:]

    return data, ref_data


def _read_hdr(filename):
    '''Read header component of a p-file.'''
    with open(filename, 'rb') as f:

        # RDB version
        f.seek(0)
        rdb_rev = np.fromfile(f, dtype=np.float32, count=1).item()
        byte_off, interp_off = _rev_info(rdb_rev)

        hdr_vars = {}
        f.seek(byte_off['off_data'])
        hdr_vars['off_data'] = np.fromfile(f, dtype=np.int32, count=1).item()

        f.seek(byte_off['nechoes'])
        hdr_vars['nechoes'] = np.fromfile(f, dtype=np.int16, count=1).item()

        f.seek(byte_off['nframes'])
        hdr_vars['nframes'] = np.fromfile(f, dtype=np.int16, count=1).item()

        f.seek(byte_off['frame_size'])
        hdr_vars['frame_size'] = np.fromfile(f, dtype=np.uint16, count=1).item()

        f.seek(byte_off['rcv'])
        hdr_vars['rcv'] = np.fromfile(f, dtype=np.int16, count=8)

        f.seek(byte_off['rhuser19'])
        hdr_vars['rhuser19'] = np.fromfile(f, dtype=np.float32, count=1).item()

        f.seek(byte_off['spec_width'])
        hdr_vars['spec_width'] = np.fromfile(f, dtype=np.float32, count=1).item()

        f.seek(byte_off['csi_dims'])
        hdr_vars['csi_dims'] = np.fromfile(f, dtype=np.int16, count=1).item()

        f.seek(byte_off['xcsi'])
        hdr_vars['xcsi'] = np.fromfile(f, dtype=np.int16, count=1).item()

        f.seek(byte_off['ycsi'])
        hdr_vars['ycsi'] = np.fromfile(f, dtype=np.int16, count=1).item()

        f.seek(byte_off['zcsi'])
        hdr_vars['zcsi'] = np.fromfile(f, dtype=np.int16, count=1).item()

        f.seek(byte_off['ps_mps_freq'])
        hdr_vars['ps_mps_freq'] = np.fromfile(f, dtype=np.uint32, count=1).astype(float).item()

        f.seek(0)
        count = max(interp_off['rdb_hdr_off_image'], interp_off['rdb_hdr_ps_mps_freq'])
        i_hdr_value = np.fromfile(f, dtype=np.int32, count=count)
        f.seek(i_hdr_value[interp_off['rdb_hdr_off_image'] - 1])
        t_hdr_value = np.fromfile(f, dtype=np.int32, count=count)

        o_hdr_off = i_hdr_value[interp_off['rdb_hdr_off_image'] - 1]
        f.seek(o_hdr_off)
        o_hdr_value = np.fromfile(f, dtype=np.float32, count=interp_off['brhc'] + 3)

    hdr_vars['tr'] = float(t_hdr_value[interp_off['image_tr']] / 1E6)
    hdr_vars['te'] = float(t_hdr_value[interp_off['image_te']] / 1E6)
    hdr_vars['vox_size'] = o_hdr_value[interp_off['image_user8']:interp_off['image_user8'] + 3]
    hdr_vars['vox_offset'] = o_hdr_value[interp_off['image_user11']:interp_off['image_user11'] + 3]
    hdr_vars['tlhc_RAS'] = o_hdr_value[interp_off['tlhc']:interp_off['tlhc'] + 3]
    hdr_vars['trhc_RAS'] = o_hdr_value[interp_off['trhc']:interp_off['trhc'] + 3]
    hdr_vars['brhc_RAS'] = o_hdr_value[interp_off['brhc']:interp_off['brhc'] + 3]

    return hdr_vars


def _rev_info(rev):
    '''Revision specific offsets.'''
    byte_offsets = {}
    interpreted_offsets = {}
    # This section is from spant
    if rev > 25:
        byte_offsets['hdr_rev'] = 0
        byte_offsets['off_data'] = 4
        byte_offsets['nechoes'] = 146
        byte_offsets['nframes'] = 150
        byte_offsets['frame_size'] = 156
        byte_offsets['rcv'] = 264
        byte_offsets['rhuser19'] = 356
        byte_offsets['spec_width'] = 432
        byte_offsets['csi_dims'] = 436
        byte_offsets['xcsi'] = 438
        byte_offsets['ycsi'] = 440
        byte_offsets['zcsi'] = 442
        byte_offsets['ps_mps_freq'] = 488
#         byte_offsets['te'] = 1148
    elif rev > 11 and rev < 25:
        byte_offsets['hdr_rev'] = 0
        byte_offsets['off_data'] = 1468
        byte_offsets['nechoes'] = 70
        byte_offsets['nframes'] = 74
        byte_offsets['frame_size'] = 80
        byte_offsets['rcv'] = 200
        byte_offsets['rhuser19'] = 292
        byte_offsets['spec_width'] = 368
        byte_offsets['csi_dims'] = 372
        byte_offsets['xcsi'] = 374
        byte_offsets['ycsi'] = 376
        byte_offsets['zcsi'] = 378
        byte_offsets['ps_mps_freq'] = 424
#         byte_offsets['te'] = 1212
    else:
        raise ValueError(f'Revision number {rev} not recognised.')

    # This section is from gannet and coroborated from pfile-tools.
    if rev == 14.3:
        interpreted_offsets['rdb_hdr_off_image']   = 377
        interpreted_offsets['rdb_hdr_ps_mps_freq'] = 107
        interpreted_offsets['image_user8']         = 37
        interpreted_offsets['image_user11']        = 40
        interpreted_offsets['tlhc']                = 120
        interpreted_offsets['trhc']                = 123
        interpreted_offsets['brhc']                = 126
        interpreted_offsets['image_te']            = 180
        interpreted_offsets['image_tr']            = 178
    elif rev == 16.0:
        interpreted_offsets['rdb_hdr_off_image']   = 377
        interpreted_offsets['rdb_hdr_ps_mps_freq'] = 107
        interpreted_offsets['image_user8']         = 49
        interpreted_offsets['image_user11']        = 52
        interpreted_offsets['tlhc']                = 132
        interpreted_offsets['trhc']                = 135
        interpreted_offsets['brhc']                = 138
        interpreted_offsets['image_te']            = 192
        interpreted_offsets['image_tr']            = 190
    elif rev in [20.006, 20.007, 24.0]:
        interpreted_offsets['rdb_hdr_off_image']   = 377
        interpreted_offsets['rdb_hdr_ps_mps_freq'] = 107
        interpreted_offsets['image_user8']         = 97
        interpreted_offsets['image_user11']        = 100
        interpreted_offsets['tlhc']                = 180
        interpreted_offsets['trhc']                = 183
        interpreted_offsets['brhc']                = 186
        interpreted_offsets['image_te']            = 266
        interpreted_offsets['image_tr']            = 264
    elif rev == 26.002:
        interpreted_offsets['rdb_hdr_off_image']   = 11
        interpreted_offsets['rdb_hdr_ps_mps_freq'] = 123
        interpreted_offsets['image_user8']         = 97
        interpreted_offsets['image_user11']        = 100
        interpreted_offsets['tlhc']                = 180
        interpreted_offsets['trhc']                = 183
        interpreted_offsets['brhc']                = 186
        interpreted_offsets['image_te']            = 266
        interpreted_offsets['image_tr']            = 264
    else:
        raise ValueError(f'Revision number {rev} not recognised.')

    return byte_offsets, interpreted_offsets
