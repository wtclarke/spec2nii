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
svs_path_sdat = philips_path / 'P1' / 'SV_PRESS_sh_6_2_raw_act.SDAT'
svs_path_spar = philips_path / 'P1' / 'SV_PRESS_sh_6_2_raw_act.SPAR'


def test_svs(tmp_path):

    subprocess.check_call(['spec2nii', 'philips',
                           '-f', 'svs',
                           '-o', tmp_path,
                           '-j',
                           str(svs_path_sdat),
                           str(svs_path_spar)])

    img_t = read_nifti_mrs(tmp_path / 'svs.nii.gz')

    assert img_t.shape == (1, 1, 1, 2048)
    assert np.iscomplexobj(img_t.dataobj)
    assert 1 / img_t.header['pixdim'][4] == 5000.0

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['SpectrometerFrequency'][0] == 127.759464
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == svs_path_sdat.name
