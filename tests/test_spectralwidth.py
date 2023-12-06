'''Simple test for SpectralWidth compliance'''

import subprocess
from pathlib import Path
import json

import numpy as np

from .io_for_tests import read_nifti_mrs

# Data paths
siemens_path = Path(__file__).parent / 'spec2nii_test_data' / 'Siemens'
vb_path = siemens_path / 'VBData' / 'Twix/meas_MID151_svs_se_C_T15_S10_10_FID108741.dat'


def test_insertion_spectralwidth(tmp_path):
    subprocess.check_call(['spec2nii', 'twix',
                           '-e', 'image',
                           '-f', 'vb',
                           '-o', tmp_path,
                           '-j', str(vb_path)])

    img_t = read_nifti_mrs(tmp_path / 'vb.nii.gz')

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert 'SpectralWidth' in hdr_ext
    assert np.isclose(hdr_ext['SpectralWidth'], 1 / img_t.header['pixdim'][4])
