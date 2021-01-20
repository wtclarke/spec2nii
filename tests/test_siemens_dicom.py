'''Tests for Siemens dicom format conversion.

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
vb_svs_path = siemens_path / 'VBData' / 'DICOM/svs_se_C>T15>S10_10_12_1/'
vb_mrsi_path = siemens_path / 'VBData' / 'DICOM/csi_se_3D_C>S23.5>T20.3_10_8_1/'

ve_svs_path = siemens_path / 'VEData' / 'DICOM/svs_se_s>c10>t5_R10_11_1/'
ve_mrsi_path = siemens_path / 'VEData' / 'DICOM/csi_se_3D_t>c23.5>s20.3_R10_7_1/'


def test_VB_svs(tmp_path):

    subprocess.check_call(['spec2nii', 'dicom',
                           '-f', 'vb_svs',
                           '-o', tmp_path,
                           '-j', str(vb_svs_path)])

    img_t = read_nifti_mrs(tmp_path / 'vb_svs.nii.gz')

    # hdr_ext_codes = img_t.header.extensions.get_codes()
    # hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_t.shape == (1, 1, 1, 512)
    assert np.iscomplexobj(img_t.dataobj)


def test_VE_svs(tmp_path):

    subprocess.check_call(['spec2nii', 'dicom',
                           '-f', 've_svs',
                           '-o', tmp_path,
                           '-j', str(ve_svs_path)])

    img_t = read_nifti_mrs(tmp_path / 've_svs.nii.gz')

    # hdr_ext_codes = img_t.header.extensions.get_codes()
    # hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_t.shape == (1, 1, 1, 512)
    assert np.iscomplexobj(img_t.dataobj)


def test_VB_mrsi(tmp_path):

    subprocess.check_call(['spec2nii', 'dicom',
                           '-f', 'vb_mrsi',
                           '-o', tmp_path,
                           '-j', str(vb_mrsi_path)])

    img_t = read_nifti_mrs(tmp_path / 'vb_mrsi.nii.gz')

    # hdr_ext_codes = img_t.header.extensions.get_codes()
    # hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_t.shape == (16, 16, 8, 512)
    assert np.iscomplexobj(img_t.dataobj)


def test_VE_mrsi(tmp_path):

    subprocess.check_call(['spec2nii', 'dicom',
                           '-f', 've_mrsi',
                           '-o', tmp_path,
                           '-j', str(ve_mrsi_path)])

    img_t = read_nifti_mrs(tmp_path / 've_mrsi.nii.gz')

    # hdr_ext_codes = img_t.header.extensions.get_codes()
    # hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_t.shape == (16, 16, 16, 512)
    assert np.iscomplexobj(img_t.dataobj)
