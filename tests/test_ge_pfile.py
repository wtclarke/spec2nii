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
mrsi_path = ge_path / 'pFiles' / 'MRSI' / 'P18432.7'
mpress_path = ge_path / 'pFiles' / 'big_gaba' / 'S01_GABA_68.7'


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

    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_DYN'
    assert hdr_ext['SpectrometerFrequency'][0] == 127.76365
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == svs_path.name

    assert np.isclose(hdr_ext['EchoTime'], 0.27)
    assert np.isclose(hdr_ext['RepetitionTime'], 2.0)


def test_mrsi(tmp_path):

    subprocess.check_call(['spec2nii', 'ge',
                           '-f', 'mrsi',
                           '-o', tmp_path,
                           '-j',
                           str(mrsi_path)])

    img = read_nifti_mrs(tmp_path / 'mrsi.nii.gz')

    assert img.shape == (8, 8, 1, 1024, 32)
    assert np.iscomplexobj(img.dataobj)
    assert 1 / img.header['pixdim'][4] == 4000.0

    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['SpectrometerFrequency'][0] == 127.763607
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == mrsi_path.name


def test_mpress(tmp_path):

    subprocess.check_call(['spec2nii', 'ge',
                           '-f', 'mpress',
                           '-o', tmp_path,
                           '-j',
                           str(mpress_path)])

    img = read_nifti_mrs(tmp_path / 'mpress.nii.gz')
    img_ref = read_nifti_mrs(tmp_path / 'mpress_ref.nii.gz')

    assert img.shape == (1, 1, 1, 2048, 8, 160, 2)
    assert np.iscomplexobj(img.dataobj)
    assert 1 / img.header['pixdim'][4] == 2000.0

    assert img_ref.shape == (1, 1, 1, 2048, 8, 8, 2)
    assert np.iscomplexobj(img_ref.dataobj)
    assert 1 / img_ref.header['pixdim'][4] == 2000.0

    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_DYN'
    assert hdr_ext['dim_7'] == 'DIM_EDIT'
    assert hdr_ext['SpectrometerFrequency'][0] == 127.758139
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == mpress_path.name

    assert np.isclose(hdr_ext['EchoTime'], 0.068)
    assert np.isclose(hdr_ext['RepetitionTime'], 2.0)
