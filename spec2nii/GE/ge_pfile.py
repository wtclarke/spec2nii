'''spec2nii module containing functions specific to interpreting the GE pfile format

Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford

This code is adapted from the VESPA project https://scion.duhs.duke.edu/vespa/project.
I therefore include their BSD statement here.

    Copyright (c) 2010, Duke University. All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
        * Redistributions of source code must retain the above copyright notice,
          this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
          notice, this list of conditions and the following disclaimer in the
          documentation and/or other materials provided with the distribution.
        * Neither the name of Duke University nor the United States Department
          of Veterans Affairs may be used to endorse or promote products derived
          from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ''AS
    IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
    THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
    LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
    OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
    DAMAGE.

    For portions of this code, copyright and license information differs from
    the above. In these cases, copyright and/or license information is inline.
'''

# Python modules
from datetime import datetime
from os.path import basename
from warnings import warn

# 3rd party modules
import numpy as np

from nifti_mrs.create_nmrs import gen_nifti_mrs_hdr_ext
from nifti_mrs.hdr_ext import Hdr_Ext

# Our modules
from spec2nii.GE.ge_read_pfile import Pfile
from spec2nii.nifti_orientation import NIFTIOrient
from spec2nii import __version__ as spec2nii_ver


class UnsupportedPulseSequenceError(Exception):
    """Raised when a user selected file does not contain the desired pulse sequence """
    pass


def read_pfile(filename, fname_out):
    """
    Given the name of a GE Pfile (*.7) file return a reference and metabolite
    data object.

    The mapper() method returns data as (x, y, z, navg, ncoil, nspectralpts)
    """

    pfile = Pfile(filename)

    if pfile.is_svs:
        data, fname_suffix = _process_svs_pfile(pfile)
    else:
        data, fname_suffix = _process_mrsi_pfile(pfile)

    fnames = []
    for fns in fname_suffix:
        if fname_out is None:
            fnames.append(filename.stem + fns)
        else:
            fnames.append(fname_out + fns)

    return data, fnames


def _process_svs_pfile(pfile):
    '''Handle SVS data

    :param Pfile pfile: Pfile object
    :return: List of NIFTI MRS data objects
    :return: List of file name suffixes
    '''
    psd = pfile.hdr.rhi_psdname.decode('utf-8').lower()

    if psd in ('probe-p', 'probe-s'):
        data, meta, dwelltime, fname_suffix = _process_probe_p(pfile)
    elif psd in ('oslaser', 'slaser_cni'):
        data, meta, dwelltime, fname_suffix = _process_oslaser(pfile)
    elif psd == 'gaba':
        data, meta, dwelltime, fname_suffix = _process_gaba(pfile)
    elif 'jpress_ac' in psd:  # Bergen patch
        data, meta, dwelltime, fname_suffix = _process_gaba(pfile)
    else:
        raise UnsupportedPulseSequenceError(f'Unrecognised sequence {psd}.')

    orientation = NIFTIOrient(_calculate_affine(pfile))

    out_nmrs = []
    for dd, mm in zip(data, meta):
        out_nmrs.append(gen_nifti_mrs_hdr_ext(dd, dwelltime, mm, orientation.Q44, no_conj=True))

    return out_nmrs, fname_suffix


def _process_probe_p(pfile):
    '''Extract metabolite and reference data from a prob_p format pfile

    :param Pfile pfile: Pfile object
    :return: List numpy data arrays
    :return: List of file name suffixes
    '''

    metab = pfile.map.raw_suppressed                    # typically (1,1,1,navg,ncoil,npts)
    metab = np.transpose(metab, [0, 1, 2, 5, 4, 3])     # swap to (1,1,1,npts,ncoil,navg)

    water = pfile.map.raw_unsuppressed                  # typically (1,1,1,navg,ncoil,npts)
    water = np.transpose(water, [0, 1, 2, 5, 4, 3])     # swap to (1,1,1,npts,ncoil,navg)

    dwelltime = 1 / pfile.hdr.rhr_spectral_width

    meta = _populate_metadata(pfile, water_suppressed=True)
    meta_ref = _populate_metadata(pfile, water_suppressed=False)

    return [metab, water], [meta, meta_ref], dwelltime, ['', '_ref']


def _process_oslaser(pfile):
    '''Extract metabolite and reference data from a oslaser format pfile

    I think this is like the CMRR sLASER sequence with 2 ecc and 2 water scaling
    scans at the start and end of each acquisition.

    :param Pfile pfile: Pfile object
    :return: List numpy data arrays
    :return: List of file name suffixes
    '''

    water = pfile.map.raw_unsuppressed                  # typically (1,1,1,navg,ncoil,npts)
    metab = pfile.map.raw_suppressed
    metab = np.transpose(metab, [0, 1, 2, 5, 4, 3])     # swap to (1,1,1,npts,ncoil,navg)
    water = np.transpose(water, [0, 1, 2, 5, 4, 3])     # swap to (1,1,1,npts,ncoil,navg)

    data = []
    meta = []
    fname_suffix = []

    data.append(water[:, :, :, :, :, [0, 1, 4, 5]])
    meta.append(_populate_metadata(pfile, water_suppressed=False))
    fname_suffix.append('_quant')
    data.append(water[:, :, :, :, :, [2, 3, 6, 7]])
    meta.append(_populate_metadata(pfile, water_suppressed=False))
    fname_suffix.append('_ecc')

    data.append(metab)
    meta.append(_populate_metadata(pfile, water_suppressed=True))
    fname_suffix.append('')

    dwelltime = 1 / pfile.hdr.rhr_spectral_width

    return data, meta, dwelltime, fname_suffix


def _process_gaba(pfile):
    '''Extract metabolite and reference data from a gaba (MPRESS) format pfile

    :param Pfile pfile: Pfile object
    :return: List numpy data arrays
    :return: List of file name suffixes
    '''

    # Note that custom mapper sorts dimensions already
    metab = pfile.map.raw_suppressed
    water = pfile.map.raw_unsuppressed

    dwelltime = 1 / pfile.hdr.rhr_spectral_width

    meta = _populate_metadata(pfile, water_suppressed=True)
    meta_ref = _populate_metadata(pfile, water_suppressed=False)

    meta.set_dim_info(0, 'DIM_COIL')
    meta.set_dim_info(1, 'DIM_DYN')
    meta.set_dim_info(2, 'DIM_EDIT')

    meta_ref.set_dim_info(0, 'DIM_COIL')
    meta_ref.set_dim_info(1, 'DIM_DYN')
    meta_ref.set_dim_info(2, 'DIM_EDIT')

    return [metab, water], [meta, meta_ref], dwelltime, ['', '_ref']


def _process_mrsi_pfile(pfile):
    '''Handle MRSI data

    :param Pfile pfile: Pfile object
    :return: List of NIFTI MRS data objects
    :return: List of file name suffixes
    '''
    psd = pfile.hdr.rhi_psdname.decode('utf-8').lower()

    known_formats = ('probe-p', 'probe-sl', 'slaser_cni', 'presscsi')
    if psd not in known_formats:
        raise UnsupportedPulseSequenceError(
            f"Unrecognised sequence {psd}, psdname must be in: {','.join(known_formats)}.")

    warn('The interpretation of pfile CSI data is poorly tested; rotations or transpositions of the'
         ' CSI grid could be present. Spec2nii currently lacks a complete set of CSI test data.'
         ' Please get in touch to help solve this issue!')

    data = np.transpose(pfile.map.raw_data, [0, 1, 2, 5, 4, 3])
    if data.shape[5] == 1:
        data = data.squeeze(axis=5)

    # Perform fft
    def fft_and_shift(x, axis):
        return np.fft.fftshift(np.fft.fft(x, axis=axis), axes=axis)

    data = fft_and_shift(data, 0)
    data = fft_and_shift(data, 1)
    data = fft_and_shift(data, 2)

    dwelltime = 1 / pfile.hdr.rhr_spectral_width
    meta = _populate_metadata(pfile)
    orientation = NIFTIOrient(_calculate_affine_mrsi(pfile))

    return [gen_nifti_mrs_hdr_ext(data, dwelltime, meta, orientation.Q44, no_conj=True), ], ['', ]


def _calculate_affine_mrsi(pfile):
    '''Calculate the 4x4 affine matrix for mrsi'''

    dcos = pfile.map.get_dcos.T
    dcos[dcos == 0.0] = 0.0             # remove -0.0 values

    voxel_size = pfile.map.get_voxel_spacing
    voxel_size_eye = np.diag(voxel_size)
    voi_position = pfile.map.get_origin_from_center(pfile.map.get_select_box_center,
                                                    pfile.map.get_num_voxels,
                                                    pfile.map.get_voxel_spacing,
                                                    pfile.map.get_dcos)

    affine = np.zeros((4, 4), dtype=float)
    affine[:3, :3] = dcos @ voxel_size_eye

    affine[:3, 3] = voi_position
    affine[3, 3] = 1.0

    return affine


def _calculate_affine(pfile):
    '''Calculate the 4x4 affine matrix'''

    dcos = pfile.map.get_dcos.T
    dcos[dcos == 0.0] = 0.0             # remove -0.0 values

    voxel_size = pfile.map.get_select_box_size
    voxel_size_eye = np.diag(voxel_size)
    voi_position = pfile.map.get_select_box_center

    affine = np.zeros((4, 4), dtype=float)
    affine[:3, :3] = dcos @ voxel_size_eye

    affine[:3, 3] = voi_position
    affine[3, 3] = 1.0

    return affine


def _populate_metadata(pfile, water_suppressed=True):
    ''' Populate a nifti-mrs header extension with the requisite information'''
    hdr = pfile.hdr
    spec_frequency = float(pfile.hdr.rhr_rh_ps_mps_freq) / 1e7

    #  Use the mps freq and field strength to determine gamma and thus isotope
    try:
        gamma = (hdr.rhr_rh_ps_mps_freq * 1e-7) / (hdr.rhe_magstrength / 10000.0)
        if abs(gamma - 42.57) < 0.3:
            nucleus = "1H"
        elif abs(gamma - 10.7) < 0.3:
            nucleus = "13C"
        elif abs(gamma - 17.2) < 0.3:
            nucleus = "31P"
        else:
            print('Warning: Unrecognised nucleus, setting to 1H as default. '
                  'Use the --override_nucleus option to specify a different nuclide.')
            nucleus = "1H"
    except ZeroDivisionError:
        # Catch data (anonymised?) which has had rhe_magstrength set to 0
        # E.g. the BIG GABA data.
        print('Warning: Magnetic field value set to zero, setting nucleus to 1H as default. '
              'Use the --override_nucleus option to specify a different nuclide.')
        nucleus = "1H"

    meta = Hdr_Ext(
        spec_frequency,
        nucleus)

    # Standard defined metadata
    # # 5.1 MRS specific Tags
    # 'EchoTime'
    meta.set_standard_def('EchoTime', float(hdr.rhi_te) / 1E6)
    # 'RepetitionTime'
    meta.set_standard_def('RepetitionTime', float(hdr.rhi_tr) / 1E6)
    # 'InversionTime'
    meta.set_standard_def('InversionTime', float(hdr.rhi_ti) / 1E6)
    # 'MixingTime'
    # Not known
    # 'ExcitationFlipAngle'
    meta.set_standard_def('ExcitationFlipAngle', float(hdr.rhi_mr_flip))
    # 'TxOffset'
    # Not known
    # 'VOI'
    # Not known
    # 'WaterSuppressed'
    meta.set_standard_def('WaterSuppressed', water_suppressed)
    # 'WaterSuppressionType'
    # Not Known
    # 'SequenceTriggered'
    # Not Known

    # # 5.2 Scanner information
    # 'Manufacturer'
    meta.set_standard_def('Manufacturer', 'GE')
    # 'ManufacturersModelName'
    meta.set_standard_def('ManufacturersModelName', hdr.rhe_ex_sysid.decode('utf-8'))
    # 'DeviceSerialNumber'
    meta.set_standard_def('DeviceSerialNumber', hdr.rhe_uniq_sys_id.decode('utf-8'))
    # 'SoftwareVersions'
    meta.set_standard_def('SoftwareVersions', hdr.rhe_ex_verscre.decode('utf-8'))
    # 'InstitutionName'
    meta.set_standard_def('InstitutionName', hdr.rhe_hospname.decode('utf-8'))
    # 'InstitutionAddress'
    # Not known
    # 'TxCoil'
    # Not Known
    # 'RxCoil'
    meta.set_user_def(key='ReceiveCoilName', value=hdr.rhi_cname.decode('utf-8'), doc='Rx coil name.')

    # # 5.3 Sequence information
    # 'SequenceName'
    meta.set_standard_def('SequenceName', hdr.rhi_psdname.decode('utf-8'))
    # 'ProtocolName'
    meta.set_standard_def('ProtocolName', hdr.rhs_se_desc.decode('utf-8'))

    # # 5.4 Sequence information
    # 'PatientPosition'
    # Not known
    # 'PatientName'
    meta.set_standard_def('PatientName', hdr.rhe_patname.decode('utf-8'))
    # 'PatientID'
    # Not known
    # 'PatientWeight'
    # Not known
    # 'PatientDoB'
    meta.set_standard_def('PatientDoB', hdr.rhe_dateofbirth.decode('utf-8'))
    # 'PatientSex'
    if hdr.rhe_patsex == 1:
        sex_str = 'M'
    elif hdr.rhe_patsex == 2:
        sex_str = 'F'
    else:
        sex_str = 'O'
    meta.set_standard_def('PatientSex', sex_str)
    # # 5.5 Provenance and conversion metadata
    # 'ConversionMethod'
    meta.set_standard_def('ConversionMethod', f'spec2nii v{spec2nii_ver}')
    # 'ConversionTime'
    conversion_time = datetime.now().isoformat(sep='T', timespec='milliseconds')
    meta.set_standard_def('ConversionTime', conversion_time)
    # 'OriginalFile'
    meta.set_standard_def('OriginalFile', [basename(pfile.file_name)])
    # # 5.6 Spatial information
    # 'kSpace'
    meta.set_standard_def('kSpace', [False, False, False])

    return meta
