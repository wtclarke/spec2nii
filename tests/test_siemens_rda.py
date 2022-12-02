'''Tests for Siemens rda format conversion.

Copyright William Clarke, University of Oxford 2022
Subject to the BSD 3-Clause License.
'''

import subprocess
from pathlib import Path
import json

import numpy as np

from .io_for_tests import read_nifti_mrs

# Data paths
siemens_path = Path(__file__).parent / 'spec2nii_test_data' / 'Siemens'

xa20_svs_path = siemens_path / 'XAData' / 'XA20/rda/spct_002.MR.MRI-LAB Test_Dir.5.1.114540.rda'
xa31_locale_svs_path = siemens_path / 'XAData' / 'XA31/rda/locale_XA31.rda'


def test_xa20_svs(tmp_path):

    subprocess.check_call(['spec2nii', 'rda',
                           '-f', 'xa_svs',
                           '-o', tmp_path,
                           '-j', str(xa20_svs_path)])

    img_t = read_nifti_mrs(tmp_path / 'xa_svs.nii.gz')

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    hdr_ext['ResonantNucleus'] = ['1H', ]

    assert img_t.shape == (1, 1, 1, 1024)
    assert np.iscomplexobj(img_t.dataobj)


def test_xa31_locale_svs(tmp_path):

    subprocess.check_call(['spec2nii', 'rda',
                           '-f', 'xa_svs',
                           '-o', tmp_path,
                           '-j', str(xa31_locale_svs_path)])

    img_t = read_nifti_mrs(tmp_path / 'xa_svs.nii.gz')

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    hdr_ext['ResonantNucleus'] = ['1H', ]

    assert img_t.shape == (1, 1, 1, 1024)
    assert np.iscomplexobj(img_t.dataobj)
