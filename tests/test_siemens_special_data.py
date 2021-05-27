'''Test special case siemens data.
Includes:
1. Twix SPECIAL (rm_special) data that doesn't report a
    twixObj.hdr.Meas.lFinalMatrixSizeSlice field.
2. Anonymised Siemens DICOM data lacking a number of normal tags.

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
data_path_special = siemens_path / 'special' / 'sub-brainphan_nuc-1H_loc-phan_spec_rm-special.dat'
data_path_anon = siemens_path / 'anon' / 'anon_dcm.IMA'


def test_twix_special(tmp_path):

    subprocess.check_call(['spec2nii', 'twix',
                           '-e', 'image',
                           '-f', 'special',
                           '-o', tmp_path,
                           '-t6', 'DIM_USER_0',
                           '-j', str(data_path_special)])

    img_t = read_nifti_mrs(tmp_path / 'special.nii.gz')

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_t.shape == (1, 1, 1, 4096, 32, 32)
    assert np.iscomplexobj(img_t.dataobj)

    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_USER_0'


def test_dicom_anon(tmp_path):

    subprocess.check_call(['spec2nii', 'dicom',
                           '-f', 'anon',
                           '-o', tmp_path,
                           '-j', str(data_path_anon)])

    img_t = read_nifti_mrs(tmp_path / 'anon.nii.gz')

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_t.shape == (1, 1, 1, 1024)
    assert np.iscomplexobj(img_t.dataobj)

    assert hdr_ext['SpectrometerFrequency'][0] == 123.255582
    assert hdr_ext['ResonantNucleus'][0] == "1H"
    assert hdr_ext['OriginalFile'][0] == str(data_path_anon)
