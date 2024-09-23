"""spec2nii module containing functions specific to interpreting the GE pfile format

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
"""

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
    """Handle SVS data

    :param Pfile pfile: Pfile object
    :return: List of NIFTI MRS data objects
    :return: List of file name suffixes
    """
    psd = pfile.hdr.rhi_psdname.decode('utf-8').lower()
    proto = pfile.hdr.rhs_se_desc.decode('utf-8').lower()
    if psd == 'hbcd' and "press" in proto:
        print('\nPSD was: ', psd)
        print('Proto is: ', proto)
        psd = pfile.hdr.rhs_se_desc.decode('utf-8').lower()
        print('PSD updated to: ', psd)

    # MM: Some 'gaba' psd strings contain full path names, so truncate to the end of the path
    if psd.endswith('gaba'):
        psd = 'gaba'

    numecho = pfile.hdr.rhi_numecho

    if psd in ('mrs-press', 'probe-p', 'probe-s', 'probe-p_ach'):
        data, meta, dwelltime, fname_suffix = _process_probe_p(pfile)
    elif psd in ('oslaser', 'slaser_cni') and numecho == 1:  # MM: If non-edited data, use _process_oslaser
        data, meta, dwelltime, fname_suffix = _process_oslaser(pfile)
    elif psd == 'oslaser' and numecho > 1:  # MM: If edited data, use _process_gaba
        data, meta, dwelltime, fname_suffix = _process_gaba(pfile)
    elif psd == 'slaser':
        data, meta, dwelltime, fname_suffix = _process_slaser(pfile)
    elif psd in ('hbcd', 'jpress', 'jpress_ac', 'gaba', 'probe-p-mega_rml', 'repress7'):
        data, meta, dwelltime, fname_suffix = _process_gaba(pfile)
    elif psd in ('hbcd2'):                                     					# ATG
        data, meta, dwelltime, fname_suffix = _process_hbcd(pfile)              # ATG
    else:
        raise UnsupportedPulseSequenceError(f'Unrecognised sequence {psd}.')

    orientation = NIFTIOrient(_calculate_affine(pfile))

    out_nmrs = []
    for dd, mm in zip(data, meta):
        out_nmrs.append(gen_nifti_mrs_hdr_ext(dd, dwelltime, mm, orientation.Q44, no_conj=True))

    return out_nmrs, fname_suffix


def _process_probe_p(pfile):
    """Extract metabolite and reference data from a prob_p format pfile

    :param Pfile pfile: Pfile object
    :return: List numpy data arrays
    :return: List of file name suffixes
    """

    metab = pfile.map.raw_suppressed                    # typically (1,1,1,navg,ncoil,npts)
    metab = np.transpose(metab, [0, 1, 2, 5, 4, 3])     # swap to (1,1,1,npts,ncoil,navg)

    water = pfile.map.raw_unsuppressed                  # typically (1,1,1,navg,ncoil,npts)
    water = np.transpose(water, [0, 1, 2, 5, 4, 3])     # swap to (1,1,1,npts,ncoil,navg)

    dwelltime = 1 / pfile.hdr.rhr_spectral_width

    meta = _populate_metadata(pfile, water_suppressed=True)
    meta_ref = _populate_metadata(pfile, water_suppressed=False)

    meta.set_dim_info(0, 'DIM_COIL')
    meta.set_dim_info(1, 'DIM_DYN')

    meta_ref.set_dim_info(0, 'DIM_COIL')
    meta_ref.set_dim_info(1, 'DIM_DYN')

    return [metab, water], [meta, meta_ref], dwelltime, ['', '_ref']


def _process_oslaser(pfile):
    """Extract metabolite and reference data from a oslaser format pfile

    I think this is like the CMRR sLASER sequence with 2 ecc and 2 water scaling
    scans at the start and end of each acquisition.

    :param Pfile pfile: Pfile object
    :return: List numpy data arrays
    :return: List of file name suffixes
    """

    water = pfile.map.raw_unsuppressed                  # typically (1,1,1,navg,ncoil,npts)
    metab = pfile.map.raw_suppressed
    metab = np.transpose(metab, [0, 1, 2, 5, 4, 3])     # swap to (1,1,1,npts,ncoil,navg)
    water = np.transpose(water, [0, 1, 2, 5, 4, 3])     # swap to (1,1,1,npts,ncoil,navg)

    data = []
    meta = []
    fname_suffix = []
    data.append(water[:, :, :, :, :, [0, 1, 4, 5]])
    meta.append(_populate_metadata(pfile, water_suppressed=False, data_dimensions=data[0].ndim))
    fname_suffix.append('_quant')
    data.append(water[:, :, :, :, :, [2, 3, 6, 7]])
    meta.append(_populate_metadata(pfile, water_suppressed=False, data_dimensions=data[1].ndim))
    fname_suffix.append('_ecc')

    data.append(metab)
    meta.append(_populate_metadata(pfile, water_suppressed=True, data_dimensions=metab.ndim))
    fname_suffix.append('')

    dwelltime = 1 / pfile.hdr.rhr_spectral_width

    return data, meta, dwelltime, fname_suffix


def _process_slaser(pfile):
    """Extract metabolite and reference data from a slaser format pfile

    This seems to be like a standard probe-p. Maybe slaser is the canonical vendor implementation.

    :param Pfile pfile: Pfile object
    :return: List numpy data arrays
    :return: List of file name suffixes
    """

    metab = pfile.map.raw_suppressed                    # typically (1,1,1,navg,ncoil,npts)
    metab = np.transpose(metab, [0, 1, 2, 5, 4, 3])     # swap to (1,1,1,npts,ncoil,navg)

    water = pfile.map.raw_unsuppressed                  # typically (1,1,1,navg,ncoil,npts)
    water = np.transpose(water, [0, 1, 2, 5, 4, 3])     # swap to (1,1,1,npts,ncoil,navg)

    dwelltime = 1 / pfile.hdr.rhr_spectral_width

    meta = _populate_metadata(pfile, water_suppressed=True, data_dimensions=metab.ndim)
    meta_ref = _populate_metadata(pfile, water_suppressed=False, data_dimensions=water.ndim)

    return [metab, water], [meta, meta_ref], dwelltime, ['', '_ref']


def _add_editing_info(pfile, meta, data):
    """Add editing information to dimension tags and headers

    :param pfile: p-file object
    :type pfile: Pfile
    :param meta: Header extension object
    :type meta: Hdr_Ext
    :param data: Shaped complex data
    :type data: np.ndarray
    """
    edit_rf_waveform = pfile.hdr.rhi_user19
    # edit_rf_waveform == 19.0 is used by HERMES and HERCULES
    if data.shape[-1] == 2 and not edit_rf_waveform == 19.0:
        edit_rf_freq_off1 = pfile.hdr.rhi_user20
        edit_rf_freq_off2 = pfile.hdr.rhi_user21
        edit_rf_ppm_off1 = edit_rf_freq_off1 / float(pfile.hdr.rhr_rh_ps_mps_freq * 1E-7)
        edit_rf_ppm_off2 = edit_rf_freq_off2 / float(pfile.hdr.rhr_rh_ps_mps_freq * 1E-7)
        edit_rf_dur = pfile.hdr.rhi_user22
        # check for default value (-1) of pulse length
        if edit_rf_dur <= 0:
            edit_rf_dur = 16000
        dim_info = "MEGA-EDITED j-difference editing, two conditions"
        dim_header = {"EditCondition": ["ON", "OFF"]}
        edit_pulse_val = {
            "ON": {"PulseOffset": edit_rf_ppm_off1, "PulseDuration": edit_rf_dur / 1E6},
            "OFF": {"PulseOffset": edit_rf_ppm_off2, "PulseDuration": edit_rf_dur / 1E6}}

        meta.set_dim_info(2, 'DIM_EDIT', hdr=dim_header, info=dim_info)
        meta.set_standard_def("EditPulse", edit_pulse_val)
    else:
        meta.set_dim_info(2, 'DIM_EDIT')


def _process_gaba(pfile):
    """Extract metabolite and reference data from a gaba (MPRESS) format pfile

    :param Pfile pfile: Pfile object
    :return: List numpy data arrays
    :return: List of file name suffixes
    """

    # Note that custom mapper sorts dimensions already
    metab = pfile.map.raw_suppressed
    water = pfile.map.raw_unsuppressed

    dwelltime = 1 / pfile.hdr.rhr_spectral_width

    meta = _populate_metadata(pfile, water_suppressed=True)
    meta_ref = _populate_metadata(pfile, water_suppressed=False)

    meta.set_dim_info(0, 'DIM_COIL')
    meta.set_dim_info(1, 'DIM_DYN')
    # Only set an EDIT dim if there is an editing dimension
    if metab.ndim == 7:
        _add_editing_info(pfile, meta, metab)

    meta_ref.set_dim_info(0, 'DIM_COIL')
    meta_ref.set_dim_info(1, 'DIM_DYN')
    # Only set an EDIT dim if there is an editing dimension
    if water.ndim == 7:
        _add_editing_info(pfile, meta_ref, water)

    return [metab, water], [meta, meta_ref], dwelltime, ['', '_ref']


def _process_hbcd(pfile):
    """
    Input:
        Pfile Object

    Output:
        List of NumPy Data Arrays
        List of File Name Suffixes

    Details:
        Hyper/ISTHMUS Sequence

        The Integrated Short-TE and Hadamard-edited Multi-Sequence (ISTHMUS)
          incorporates a Short TE (35ms) PRESS, Long-TE (80ms) HERCULES, and
          a water reference for each.

        Data is organized within the file as follows:
          ( 1) Long  TE Reference : 80ms Unsupressed Water
          (32) Long  TE Edited    : 80ms Water Suppressed HERCULES
          ( 1) Short TE Reference : 35ms Unsupressed Water
          (32) Short TE Unedited  : 35ms Water Suppressed PRESS
          ( 1) Long  TE Reference : 80ms Unsupressed Water
          (32) Short TE Unedited  : 35ms Water Suppressed PRESS
          ( 1) Short TE Reference : 35ms Unsupressed Water
          (32) Short TE Unedited  : 35ms Water Suppressed PRESS
          ( 1) Long  TE Reference : 80ms Unsupressed Water
          (32) Short TE Unedited  : 35ms Water Suppressed PRESS
          ( 1) Short TE Reference : 35ms Unsupressed Water
          (32) Short TE Unedited  : 35ms Water Suppressed PRESS
          ( 1) Long  TE Reference : 80ms Unsupressed Water
          (32) Short TE Unedited  : 35ms Water Suppressed PRESS
          ( 1) Short TE Reference : 35ms Unsupressed Water
          (32) Short TE Unedited  : 35ms Water Suppressed PRESS

        Data is directly separated from the raw data (pfile.map.raw_data) where the data
          mapper (GABA mapper) is simply used to populate in the raw data.

        Author : Aaron Gudmundson, Johns Hopkins University, 2024
        Contact: agudmun2@jhmi.edu
    """

    # Additional Imports
    import copy

    # Editing Parameters
    edit_cases       = 4                                                                    # 4 Editing Conditions
    edit_pulse_1     = 4.58                                                                 # 4.58 ppm
    edit_pulse_2     = 1.90                                                                 # 1.90 ppm
    edit_pulse_4     = 4.18                                                                 # 4.18 ppm
    pulse_length     = 0.02                                                                 # Edit Pulse 20 ms

    dim_header       = {'EditCondition': ['A', 'B', 'C', 'D']}                              # 4 Subscans
    edit_pulse_val   = {'A': {'PulseOffset': [edit_pulse_1, edit_pulse_2], 'PulseDuration': pulse_length},
                        'B': {'PulseOffset': [edit_pulse_4, edit_pulse_2], 'PulseDuration': pulse_length},
                        'C': {'PulseOffset': edit_pulse_1, 'PulseDuration': pulse_length},
                        'D': {'PulseOffset': edit_pulse_4, 'PulseDuration': pulse_length}}

    # All Data (Skip 1st Transient - GE automatically has historically included a 'noise' transient)
    raw_data         = pfile.map.raw_data[:, :, :, :, 1:, :]                                # Raw Data from Mapper

    # Long TE HERCULES Metabolite Data
    lTE_metab        = copy.deepcopy(raw_data)                                              # Long TE Metab
    lTE_mask         = np.ones(lTE_metab.shape[4], dtype=bool)                              # Create a Mask
    lTE_mask[::33]   = False                                                                # Remove Water Refs
    lTE_mask[: 33]   = False                                                                # Remove PRESS
    lTE_metab        = lTE_metab[:, :, :, :, lTE_mask, :]                                   # Isolated HERCULES

    # Handle Incomplete
    if lTE_mask.shape[-1] % 4 != 0:                                                         # Incomplete Acquisition
        old_num_avgs = lTE_mask.shape[-1]                                                   # Old Total Averages
        new_num_avgs = (lTE_mask.shape[-1] // 4) * 4                                        # New Total Averages
        lTE_metab    = lTE_metab[:, :, :, :, :new_num_avgs, :]                              # Remove Incomplete

        notestring   = '80ms HERCULES'                                                      # Note Incomplete Data
        notestring   = f'{notestring} - Correcting - Incomplete Averages'                   # Note Incomplete Data
        notestring   = f'{notestring}  {old_num_avgs} --> {new_num_avgs}'                   # Note Incomplete Data
        print(f'{notestring} \t Corrected**')                                               # Note Incomplete Data

    new_shape     = list(lTE_metab.shape)                                                   # Remove Averages Dim
    new_shape[4]  = new_shape[4] // 4                                                       # Closest multiple of 4
    new_shape.append(edit_cases)                                                        	# Include Subscans

    lTE_metab_    = np.zeros(new_shape, dtype=np.complex128)                                # New lTE Metab Array
    lTE_metab_[:, :, :, :, :, :, 0] = lTE_metab[:, :, :, :, 0::4, :]                        # Subscan 1
    lTE_metab_[:, :, :, :, :, :, 1] = lTE_metab[:, :, :, :, 1::4, :]                        # Subscan 2
    lTE_metab_[:, :, :, :, :, :, 2] = lTE_metab[:, :, :, :, 2::4, :]                        # Subscan 3
    lTE_metab_[:, :, :, :, :, :, 3] = lTE_metab[:, :, :, :, 3::4, :]                        # Subscan 4
    lTE_metab      = lTE_metab_                                                             # With Subscan Dim

    lTE_metab_meta = _populate_metadata(pfile, water_suppressed=True)                       # Acquisition Information
    lTE_metab_meta.set_standard_def('EchoTime', 0.080)                                      # TE
    lTE_metab_meta.set_standard_def('WaterSuppressed', True)                                # Water Suppression
    lTE_metab_meta.set_standard_def('EditPulse', edit_pulse_val)                            # Header Edit Info

    lTE_metab_meta.set_dim_info(0, 'DIM_DYN')                                               # Dimension Info
    lTE_metab_meta.set_dim_info(1, 'DIM_COIL')                                              # Dimension Info
    lTE_metab_meta.set_dim_info(2, 'DIM_EDIT', hdr=dim_header)                              # Dimension Info

    # Short TE HERCULES Metabolite Data
    sTE_metab = copy.deepcopy(raw_data[:, :, :, :, 1:33, :])

    sTE_metab_meta = _populate_metadata(pfile, water_suppressed=True)                       # Acquisition Information
    sTE_metab_meta.set_standard_def('EchoTime', 0.035)                                      # TE
    sTE_metab_meta.set_standard_def('WaterSuppressed', True)                                # Water Suppression

    sTE_metab_meta.set_dim_info(0, 'DIM_DYN')                                               # Dimension Info
    sTE_metab_meta.set_dim_info(1, 'DIM_COIL')                                              # Dimension Info

    # Long TE Reference Water Data
    lTE_water = copy.deepcopy(raw_data[:, :, :, :, 0::66, :])

    lTE_water_meta = _populate_metadata(pfile, water_suppressed=False)                      # Acquisition Information
    lTE_water_meta.set_standard_def('EchoTime', 0.080)                                      # TE
    lTE_water_meta.set_standard_def('WaterSuppressed', False)                               # Water Suppression

    lTE_water_meta.set_dim_info(0, 'DIM_DYN')                                               # Dimension Info
    lTE_water_meta.set_dim_info(1, 'DIM_COIL')                                              # Dimension Info

    # Short TE Reference Water Data
    sTE_water = copy.deepcopy(raw_data[:, :, :, :, 33::66, :])

    sTE_water_meta = _populate_metadata(pfile, water_suppressed=False)                      # Acquisition Information
    sTE_water_meta.set_standard_def('EchoTime', 0.035)                                      # TE
    sTE_water_meta.set_standard_def('WaterSuppressed', False)                               # Water Suppression

    sTE_water_meta.set_dim_info(0, 'DIM_DYN')                                               # Dimension Info
    sTE_water_meta.set_dim_info(1, 'DIM_COIL')                                              # Dimension Info

    # Dwell Time
    dwelltime = 1 / pfile.hdr.rhr_spectral_width

    data      = [lTE_metab, sTE_metab, lTE_water, sTE_water]                                # ISTHMUS Data
    meta      = [lTE_metab_meta, sTE_metab_meta, lTE_water_meta, sTE_water_meta]            # ISTHMUS Header
    ref_names = ['_edited', '_short_te', '_ref_edited', '_ref_short_te']                    # ISTHMUS Naming

    print('Returning ISTHMUS Data:')
    for ii in range(len(data)):
        print('    {:02d} {:<14} '.format(ii, ref_names[ii]), data[ii].shape)
    print(' ')

    return data, meta, dwelltime, ref_names


def _process_mrsi_pfile(pfile):
    """Handle MRSI data

    :param Pfile pfile: Pfile object
    :return: List of NIFTI MRS data objects
    :return: List of file name suffixes
    """
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
    meta.set_dim_info(0, 'DIM_COIL')

    orientation = NIFTIOrient(_calculate_affine_mrsi(pfile))

    return [gen_nifti_mrs_hdr_ext(data, dwelltime, meta, orientation.Q44, no_conj=True), ], ['', ]


def _calculate_affine_mrsi(pfile):
    """Calculate the 4x4 affine matrix for mrsi"""

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
    """Calculate the 4x4 affine matrix"""

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


def _populate_metadata(pfile, water_suppressed=True, data_dimensions=None):
    """Populate a nifti-mrs header extension with metadata from the pfile

    If (up to 7) data_dimensions are specified then default dimension tags
    (coil, dyn, indirect) will be included. Otherwise manually specify
    outside this function.

    :param pfile: pfile object
    :type pfile: pfile map object
    :param water_suppressed: Set water suppression header field, defaults to True
    :type water_suppressed: bool, optional
    :param data_dimensions: If set to 5,6, or 7 will include default dim tags for those dimensions, defaults to None
    :type data_dimensions: int, optional
    :return: Header extension object
    :rtype: nifti_mrs.hdr_ext
    """
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
        nucleus,
        dimensions=data_dimensions)

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
