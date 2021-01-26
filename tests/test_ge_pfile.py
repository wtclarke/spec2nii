'''Tests for GE pfile format conversion.

Copyright William Clarke, University of Oxford 2021
Subject to the BSD 3-Clause License.
'''

import subprocess
from pathlib import Path
import json

import numpy as np

from .io_for_tests import read_nifti_mrs

# Data paths
ge_path = Path(__file__).parent / 'spec2nii_test_data' / 'ge'
svs_path = ge_path / 'pFiles' / 'svMRS' / 'P03072.7'


def test_svs(tmp_path):

    subprocess.check_call(['spec2nii', 'ge',
                           '-f', 'svs',
                           '-o', tmp_path,
                           '-j',
                           str(svs_path)])

    img = read_nifti_mrs(tmp_path / 'svs.nii.gz')
    img_ref = read_nifti_mrs(tmp_path / 'svs_ref.nii.gz')

    assert img.shape == (1, 1, 1, 4096, 32, 2)
    assert np.iscomplexobj(img.dataobj)
    assert 1 / img.header['pixdim'][4] == 5000.0

    assert img_ref.shape == (1, 1, 1, 4096, 32, 2)
    assert np.iscomplexobj(img_ref.dataobj)
    assert 1 / img_ref.header['pixdim'][4] == 5000.0

    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['SpectrometerFrequency'][0] == 127.76365
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == svs_path.name
