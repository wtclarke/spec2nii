'''Tests for Siemens twix format conversion.

Copyright William Clarke, University of Oxford 2021
Subject to the BSD 3-Clause License.
'''

import subprocess
from pathlib import Path
import json

import numpy as np

from .io_for_tests import read_nifti_mrs

# Data paths
siemens_path = Path(__file__).parent / 'spec2nii_test_data' / 'Siemens'
vb_path = siemens_path / 'VBData' / 'Twix/meas_MID151_svs_se_C_T15_S10_10_FID108741.dat'
ve_path = siemens_path / 'VEData' / 'Twix/meas_MID00240_FID62745_svs_se_c_t15_s10_R10.dat'


def test_VB(tmp_path):

    subprocess.check_call(['spec2nii', 'twix',
                           '-e', 'image',
                           '-f', 'vb',
                           '-o', tmp_path,
                           '-j', str(vb_path)])

    img_t = read_nifti_mrs(tmp_path / 'vb.nii.gz')

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_t.shape == (1, 1, 1, 1056, 32, 2)
    assert np.iscomplexobj(img_t.dataobj)

    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_DYN'


def test_VE(tmp_path):

    subprocess.check_call(['spec2nii', 'twix',
                           '-e', 'image',
                           '-f', 've',
                           '-o', tmp_path,
                           '-j', str(ve_path)])

    img_t = read_nifti_mrs(tmp_path / 've.nii.gz')

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_t.shape == (1, 1, 1, 1056, 32, 16)
    assert np.iscomplexobj(img_t.dataobj)

    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_DYN'
