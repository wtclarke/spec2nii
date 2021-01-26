'''GE p-file reader

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

# 3rd party modules
import numpy as np

# Our modules
from spec2nii.GE.ge_read_pfile import Pfile
from spec2nii import nifti_mrs
from spec2nii.nifti_orientation import NIFTIOrient


class UnsupportedPulseSequenceError(Exception):
    """Raised when a user selected file does not contain the desired pulse sequence """
    pass


class SIDataError(Exception):
    '''Currently we don't handle SI data.'''
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
        raise SIDataError

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

    if psd == 'probe-p':
        data, fname_suffix = _process_probe_p(pfile)
    elif psd == 'oslaser':
        data, fname_suffix = _process_oslaser(pfile)
    else:
        raise UnsupportedPulseSequenceError('Unrecognised sequence, psdname must be "prob-p" or "oslaser".')

    dwelltime = 1 / pfile.hdr.rhr_spectral_width
    meta = []

    for dd in data:
        mm = _populate_metadata(pfile)
        meta.append(mm)

    orientation = NIFTIOrient(_calculate_affine(pfile))

    out_nmrs = []
    for dd, mm in zip(data, meta):
        out_nmrs.append(nifti_mrs.NIfTI_MRS(dd, orientation.Q44, dwelltime, mm))

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

    return [metab, water], ['', '_ref']


def _process_oslaser(pfile):
    '''Extract metabolite and reference data from a oslaser format pfile

    I think this is the CMRR sLASER sequence with 2 ecc and 2 water scaling
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
    fname_suffix = []

    data.append(water[:, :, :, :, :, [0, 1, 4, 5]])
    fname_suffix.append('_quant')
    data.append(water[:, :, :, :, :, [2, 3, 6, 7]])
    fname_suffix.append('_ecc')

    data.append(metab)
    fname_suffix.append('')

    return data, fname_suffix


def _calculate_affine(pfile):
    '''Calculate the 4x4 affine matrix'''

    dcos = pfile.map.get_dcos.T
    dcos[dcos == 0.0] = 0.0             # remove -0.0 values

    voxel_size = np.diag(pfile.map.get_select_box_size)
    voi_position = pfile.map.get_select_box_center

    affine = np.zeros((4, 4), dtype=np.float)
    affine[:3, :3] = dcos @ voxel_size

    affine[:3, 3] = voi_position
    affine[3, 3] = 1.0

    return affine


def _populate_metadata(pfile):
    ''' Populate a nifti-mrs header extension with the requisite information'''
    hdr = pfile.hdr
    spec_frequency = float(pfile.hdr.rhr_rh_ps_mps_freq) / 1e7

    #  Use the mps freq and field strength to determine gamma and thus isotope
    gamma = (hdr.rhr_rh_ps_mps_freq * 1e-7) / (hdr.rhe_magstrength / 10000.0)
    if abs(gamma - 42.57) < 0.3:
        nucleus = "1H"
    elif abs(gamma - 10.7) < 0.3:
        nucleus = "13C"
    elif abs(gamma - 17.2) < 0.3:
        nucleus = "31P"
    else:
        nucleus = "1H"

    meta = nifti_mrs.hdr_ext(spec_frequency,
                             nucleus)

    meta.set_standard_def('Manufacturer', 'GE')
    meta.set_standard_def('ManufacturersModelName', hdr.rhe_ex_sysid.decode('utf-8'))
    meta.set_standard_def('DeviceSerialNumber', hdr.rhe_uniq_sys_id.decode('utf-8'))
    meta.set_standard_def('SoftwareVersions', hdr.rhe_ex_verscre.decode('utf-8'))

    meta.set_standard_def('InstitutionName', hdr.rhe_hospname.decode('utf-8'))
    # meta.set_standard_def('InstitutionAddress', )

    meta.set_user_def(key='ReceiveCoilName', value=hdr.rhi_cname.decode('utf-8'), doc='Rx coil name.')

    # Some sequence information
    meta.set_standard_def('SequenceName', hdr.rhi_psdname.decode('utf-8'))
    meta.set_standard_def('ProtocolName', hdr.rhs_se_desc.decode('utf-8'))

    # Some subject information
    # meta.set_standard_def('PatientPosition',)
    meta.set_standard_def('PatientName', hdr.rhe_patname.decode('utf-8'))
    # meta.set_standard_def('PatientWeight', )
    meta.set_standard_def('PatientDoB', hdr.rhe_dateofbirth.decode('utf-8'))

    if hdr.rhe_patsex == 1:
        sex_str = 'M'
    elif hdr.rhe_patsex == 2:
        sex_str = 'F'
    else:
        sex_str = 'O'
    meta.set_standard_def('PatientSex', sex_str)

    # Timing and sequence parameters
    meta.set_standard_def('InversionTime', hdr.rhi_ti)
    meta.set_standard_def('ExcitationFlipAngle', hdr.rhi_mr_flip)
    # meta.set_standard_def('TxOffset', )
    meta.set_standard_def('EchoTime', float(hdr.rhi_te))
    meta.set_standard_def('RepetitionTime', float(hdr.rhi_tr))

    meta.set_standard_def('ConversionMethod', 'spec2nii')
    conversion_time = datetime.now().isoformat(sep='T', timespec='milliseconds')
    meta.set_standard_def('ConversionTime', conversion_time)
    meta.set_standard_def('OriginalFile', [basename(pfile.file_name)])

    return meta
