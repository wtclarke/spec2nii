'''Tests for Philips SDAT/SPAR format conversion.

Copyright William Clarke, University of Oxford 2021
Subject to the BSD 3-Clause License.
'''

import subprocess
from pathlib import Path
import json

import numpy as np

from .io_for_tests import read_nifti_mrs

# Data paths
philips_path = Path(__file__).parent / 'spec2nii_test_data' / 'philips'
list_path = philips_path / 'hyper' / 'raw_226.list'
data_path = philips_path / 'hyper' / 'raw_226.data'
aux_path = philips_path / 'hyper' / 'converted_dcm.dcm'


def test_svs_hyper(tmp_path):

    subprocess.check_call(['spec2nii', 'philips_dl',
                           '-f', 'svs',
                           '-o', tmp_path,
                           '-j',
                           '--special', 'hyper',
                           str(data_path),
                           str(list_path),
                           str(aux_path)])

    img_1 = read_nifti_mrs(tmp_path / 'svs_hyper_short_te.nii.gz')
    img_2 = read_nifti_mrs(tmp_path / 'svs_hyper_edited.nii.gz')
    img_3 = read_nifti_mrs(tmp_path / 'svs_NOI.nii.gz')

    assert img_1.shape == (1, 1, 1, 2048, 15, 32)
    assert np.iscomplexobj(img_1.dataobj)
    assert 1 / img_1.header['pixdim'][4] == 2000.0

    assert img_2.shape == (1, 1, 1, 2048, 15, 4, 56)
    assert np.iscomplexobj(img_2.dataobj)
    assert 1 / img_2.header['pixdim'][4] == 2000.0

    assert img_3.shape == (1, 1, 1, 2048, 15)
    assert np.iscomplexobj(img_2.dataobj)
    assert 1 / img_2.header['pixdim'][4] == 2000.0

    hdr_ext_codes = img_1.header.extensions.get_codes()
    hdr_ext = json.loads(img_1.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['SpectrometerFrequency'][0] == 127.74876
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == data_path.name
    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_DYN'

    hdr_ext_codes = img_2.header.extensions.get_codes()
    hdr_ext = json.loads(img_2.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_EDIT'
    assert hdr_ext['dim_6_header'] == {'EditCondition': ['A', 'B', 'C', 'D']}
    assert hdr_ext['dim_7'] == 'DIM_DYN'
