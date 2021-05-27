'''Tests for NIfTI-MRS anonymisation.

Copyright William Clarke, University of Oxford 2021
Subject to the BSD 3-Clause License.
'''

import numpy as np
import subprocess
from .io_for_tests import read_nifti_mrs
from pathlib import Path
import json

# Data paths
siemens_path = Path(__file__).parent / 'spec2nii_test_data' / 'Siemens'
data_path = siemens_path / 'VBData' / 'Twix' / 'meas_MID151_svs_se_C_T15_S10_10_FID108741.dat'


def test_anon(tmp_path):
    # Convert twix
    subprocess.check_call(['spec2nii', 'twix',
                           '-e', 'image',
                           '-f', 'original',
                           '-o', tmp_path,
                           '-j', str(data_path)])

    subprocess.check_call(['spec2nii', 'anon',
                           '-f', 'anon',
                           '-o', tmp_path,
                           '-r', 'PatientSex',
                           '-r', 'RepetitionTime',
                           '-v', str(tmp_path / 'original.nii.gz')])

    img_o = read_nifti_mrs(tmp_path / 'original.nii.gz')
    img_a = read_nifti_mrs(tmp_path / 'anon.nii.gz')

    assert img_o.shape == img_a.shape
    assert np.allclose(img_o.get_fdata(dtype=complex), img_a.get_fdata(dtype=complex))
    assert np.allclose(img_o.affine, img_a.affine)
    assert np.isclose(img_o.header['pixdim'][4], img_a.header['pixdim'][4])

    hdr_ext_codes = img_o.header.extensions.get_codes()
    hdr_ext_o = json.loads(img_o.header.extensions[hdr_ext_codes.index(44)].get_content())
    hdr_ext_codes = img_a.header.extensions.get_codes()
    hdr_ext_a = json.loads(img_a.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert 'SpectrometerFrequency' in hdr_ext_a
    assert 'ResonantNucleus' in hdr_ext_a
    assert 'dim_5' in hdr_ext_a
    assert 'OriginalFile' in hdr_ext_o
    assert 'OriginalFile' not in hdr_ext_a
    assert 'PatientName' in hdr_ext_o
    assert 'PatientName' not in hdr_ext_a
    assert 'PatientSex' not in hdr_ext_a
    assert 'RepetitionTime' not in hdr_ext_a
