'''Tests for Varian format conversion.

Copyright Jack Miller, University of Oxford and University of Aarhus, 2021
Subject to the BSD 3-Clause License.
'''

import subprocess
from pathlib import Path
import json

import numpy as np

from .io_for_tests import read_nifti_mrs

# Data paths
varian_path = Path(__file__).parent / 'spec2nii_test_data' / 'varian'
data_path = varian_path / '13C_spuls_01.fid'


def test_svs(tmp_path):

    subprocess.check_call(['spec2nii', 'varian',
                           '-f', 'svs',
                           '-o', tmp_path,
                           '-j',
                           '-d',
                           '-t', 'DIM_USER_0',
                           str(data_path)])

    img_t = read_nifti_mrs(tmp_path / 'svs.nii.gz')

    assert img_t.shape == (1, 1, 1, 2048, 1, 80)
    assert np.iscomplexobj(img_t.dataobj)
    assert np.isclose(1 / img_t.header['pixdim'][4], 40323, atol=2)

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['SpectrometerFrequency'][0] == 100.4451626
    assert hdr_ext['ResonantNucleus'][0] == "13C"
    assert hdr_ext['OriginalFile'][0] == data_path.name
    assert 'VarianProcpar' in hdr_ext
    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_USER_0'
