'''Tests for Bruker format conversion.

Copyright William Clarke, University of Oxford 2021
Subject to the BSD 3-Clause License.
'''

import subprocess
from pathlib import Path
import json

import numpy as np

from .io_for_tests import read_nifti_mrs

# Data paths
bruker_path = Path(__file__).parent / 'spec2nii_test_data' / 'bruker'
data_path = bruker_path / '20201208_105201_lego_rod_1_3'


def test_fid(tmp_path):

    subprocess.check_call(['spec2nii', 'bruker',
                           '-f', 'fid',
                           '-m', 'FID',
                           '-d',
                           '-o', tmp_path,
                           '-j',
                           str(data_path)])

    # Img 5 - csi
    img = read_nifti_mrs(tmp_path / 'fid_FID_5.nii.gz')

    assert img.shape == (16, 16, 1, 1980)
    assert np.iscomplexobj(img.dataobj)
    assert 1 / img.header['pixdim'][4] == 4000.0

    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert np.isclose(hdr_ext['SpectrometerFrequency'][0], 400.32251)
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == str(data_path.absolute() / '5' / 'fid')
    assert 'method' in hdr_ext
    assert 'acqp' in hdr_ext

    # Img 6 - csi
    img = read_nifti_mrs(tmp_path / 'fid_FID_6.nii.gz')

    assert img.shape == (16, 16, 1, 1980)
    assert np.iscomplexobj(img.dataobj)
    assert 1 / img.header['pixdim'][4] == 4000.0

    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert np.isclose(hdr_ext['SpectrometerFrequency'][0], 400.32251)
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == str(data_path.absolute() / '6' / 'fid')
    assert 'method' in hdr_ext
    assert 'acqp' in hdr_ext

    # Img 9 - svs
    img = read_nifti_mrs(tmp_path / 'fid_FID_9.nii.gz')

    assert img.shape == (1, 1, 1, 1980)
    assert np.iscomplexobj(img.dataobj)
    assert np.isclose(1 / img.header['pixdim'][4], 4401.41)

    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert np.isclose(hdr_ext['SpectrometerFrequency'][0], 400.32251)
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == str(data_path.absolute() / '9' / 'fid')
    assert 'method' in hdr_ext
    assert 'acqp' in hdr_ext

    # Img 10 - svs
    img = read_nifti_mrs(tmp_path / 'fid_FID_10.nii.gz')

    assert img.shape == (1, 1, 1, 1980)
    assert np.iscomplexobj(img.dataobj)
    assert np.isclose(1 / img.header['pixdim'][4], 4401.41)

    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert np.isclose(hdr_ext['SpectrometerFrequency'][0], 400.32251)
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == str(data_path.absolute() / '10' / 'fid')
    assert 'method' in hdr_ext
    assert 'acqp' in hdr_ext


def test_2dseq(tmp_path):

    subprocess.check_call(['spec2nii', 'bruker',
                           '-f', '2dseq',
                           '-m', '2DSEQ',
                           '-d',
                           '-o', tmp_path,
                           '-j',
                           str(data_path)])

    # Img 5 - csi
    img = read_nifti_mrs(tmp_path / '2dseq_2DSEQ_5_2_lego_rod_3.nii.gz')

    assert img.shape == (16, 16, 1, 2048)
    assert np.iscomplexobj(img.dataobj)
    assert 1 / img.header['pixdim'][4] == 4000.0

    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert np.isclose(hdr_ext['SpectrometerFrequency'][0], 400.32251)
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == str(data_path.absolute() / '5' / 'pdata' / '2' / '2dseq')
    assert 'method' in hdr_ext
    assert 'visu_pars' in hdr_ext

    # Img 6 - csi
    img = read_nifti_mrs(tmp_path / '2dseq_2DSEQ_6_2_lego_rod_3.nii.gz')

    assert img.shape == (16, 16, 1, 2048)
    assert np.iscomplexobj(img.dataobj)
    assert 1 / img.header['pixdim'][4] == 4000.0

    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert np.isclose(hdr_ext['SpectrometerFrequency'][0], 400.32251)
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == str(data_path.absolute() / '6' / 'pdata' / '2' / '2dseq')
    assert 'method' in hdr_ext
    assert 'visu_pars' in hdr_ext

    # Img 9 - svs
    img = read_nifti_mrs(tmp_path / '2dseq_2DSEQ_9_2_lego_rod_3.nii.gz')

    assert img.shape == (1, 1, 1, 2048)
    assert np.iscomplexobj(img.dataobj)
    assert np.isclose(1 / img.header['pixdim'][4], 4401.41)

    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert np.isclose(hdr_ext['SpectrometerFrequency'][0], 400.32251)
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == str(data_path.absolute() / '9' / 'pdata' / '2' / '2dseq')
    assert 'method' in hdr_ext
    assert 'visu_pars' in hdr_ext

    # Img 10 - svs
    img = read_nifti_mrs(tmp_path / '2dseq_2DSEQ_10_2_lego_rod_3.nii.gz')

    assert img.shape == (1, 1, 1, 2048)
    assert np.iscomplexobj(img.dataobj)
    assert np.isclose(1 / img.header['pixdim'][4], 4401.41)

    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert np.isclose(hdr_ext['SpectrometerFrequency'][0], 400.32251)
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == str(data_path.absolute() / '10' / 'pdata' / '2' / '2dseq')
    assert 'method' in hdr_ext
    assert 'visu_pars' in hdr_ext
