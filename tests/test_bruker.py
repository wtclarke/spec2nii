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

def test_fid_pv_6_0_1(tmp_path):
    data_path = bruker_path / 'PV-6.0.1'

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


def test_2dseq_pv_6_0_1(tmp_path):
    data_path = bruker_path / 'PV-6.0.1'

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


def test_fid_pv_360_1_1(tmp_path):
    data_path = bruker_path / 'PV-360.1.1'

    subprocess.check_call(['spec2nii', 'bruker',
                           '-f', 'fid',
                           '-m', 'FID',
                           '-d',
                           '-o', tmp_path,
                           '-j',
                           '-zp',
                           str(data_path)])

    # Img 7 - svs
    img = read_nifti_mrs(tmp_path / 'fid_FID_7.nii.gz')

    assert img.shape == (1, 1, 1, 2048)
    assert np.iscomplexobj(img.dataobj)
    assert np.isclose(1 / img.header['pixdim'][4], 3333.33)

    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert np.isclose(hdr_ext['SpectrometerFrequency'][0], 300.32)
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == str(data_path.absolute() / '7' / 'fid')
    assert 'method' in hdr_ext
    assert 'acqp' in hdr_ext

    # Img 12 is missing the fid file

    # Img 15 - svs
    img = read_nifti_mrs(tmp_path / 'fid_FID_15.nii.gz')

    assert img.shape == (1, 1, 1, 2048)
    assert np.iscomplexobj(img.dataobj)
    assert np.isclose(1 / img.header['pixdim'][4], 3333.33)

    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert np.isclose(hdr_ext['SpectrometerFrequency'][0], 300.32)
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == str(data_path.absolute() / '15' / 'fid')
    assert 'method' in hdr_ext
    assert 'acqp' in hdr_ext


def test_rawdata_pv_360_1_1(tmp_path):
    data_path = bruker_path / 'PV-360.1.1'

    subprocess.check_call(['spec2nii', 'bruker',
                           '-f', 'raw',
                           '-m', 'RAWDATA',
                           '-d',
                           '-o', tmp_path,
                           '-j',
                           '-zp',
                           str(data_path)])

    # Img 7 - svs
    img = read_nifti_mrs(tmp_path / 'raw_RAWDATA_7.nii.gz')

    assert img.shape == (1, 1, 1, 2048, 4)
    assert np.iscomplexobj(img.dataobj)
    assert np.isclose(1 / img.header['pixdim'][4], 3333.33)

    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert np.isclose(hdr_ext['SpectrometerFrequency'][0], 300.32)
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == str(data_path.absolute() / '7' / 'rawdata.job0')
    assert 'method' in hdr_ext
    assert 'acqp' in hdr_ext

    # Img 12 - svs
    img = read_nifti_mrs(tmp_path / 'raw_RAWDATA_12.nii.gz')

    assert img.shape == (1, 1, 1, 2060, 4, 384)
    assert np.iscomplexobj(img.dataobj)
    assert np.isclose(1 / img.header['pixdim'][4], 2949.85)

    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert np.isclose(hdr_ext['SpectrometerFrequency'][0], 300.32)
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == str(data_path.absolute() / '12' / 'rawdata.job0')
    assert 'method' in hdr_ext
    assert 'acqp' in hdr_ext

    # Img 15 - svs
    img = read_nifti_mrs(tmp_path / 'raw_RAWDATA_15.nii.gz')

    assert img.shape == (1, 1, 1, 2048, 4)
    assert np.iscomplexobj(img.dataobj)
    assert np.isclose(1 / img.header['pixdim'][4], 3333.33)

    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert np.isclose(hdr_ext['SpectrometerFrequency'][0], 300.32)
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == str(data_path.absolute() / '15' / 'rawdata.job0')
    assert 'method' in hdr_ext
    assert 'acqp' in hdr_ext


def test_fid_pv_360_3_5(tmp_path):
    data_path = bruker_path / 'PV-360.3.5'

    subprocess.check_call(['spec2nii', 'bruker',
                           '-f', 'fid',
                           '-m', 'FID',
                           '-d',
                           '-o', tmp_path,
                           '-j',
                           '-zp',
                           str(data_path)])

    # Img 7 - svs
    img = read_nifti_mrs(tmp_path / 'fid_FID_7_1.nii.gz')

    assert img.shape == (1, 1, 1, 2048)
    assert np.iscomplexobj(img.dataobj)
    assert np.isclose(1 / img.header['pixdim'][4], 3012.05)

    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert np.isclose(hdr_ext['SpectrometerFrequency'][0], 300.32)
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == str(data_path.absolute() / '7' / 'pdata' / '1' / 'fid_proc.64')
    assert 'method' in hdr_ext
    assert 'acqp' in hdr_ext

    # Img 19 - svs
    img = read_nifti_mrs(tmp_path / 'fid_FID_19_1.nii.gz')

    assert img.shape == (1, 1, 1, 2048, 2)
    assert np.iscomplexobj(img.dataobj)
    assert np.isclose(1 / img.header['pixdim'][4], 3012.05)

    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert np.isclose(hdr_ext['SpectrometerFrequency'][0], 300.32)
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == str(data_path.absolute() / '19' / 'pdata' / '1' / 'fid_proc.64')
    assert 'method' in hdr_ext
    assert 'acqp' in hdr_ext


def test_rawdata_pv_360_3_5(tmp_path):
    data_path = bruker_path / 'PV-360.3.5'

    subprocess.check_call(['spec2nii', 'bruker',
                           '-f', 'raw',
                           '-m', 'RAWDATA',
                           '-d',
                           '-o', tmp_path,
                           '-j',
                           '-zp',
                           str(data_path)])

    # Img 7 - svs
    img = read_nifti_mrs(tmp_path / 'raw_RAWDATA_7.nii.gz')

    assert img.shape == (1, 1, 1, 2048, 4)
    assert np.iscomplexobj(img.dataobj)
    assert np.isclose(1 / img.header['pixdim'][4], 3012.05)

    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert np.isclose(hdr_ext['SpectrometerFrequency'][0], 300.32)
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == str(data_path.absolute() / '7' / 'rawdata.job0')
    assert 'method' in hdr_ext
    assert 'acqp' in hdr_ext

    # Img 19 - svs
    img = read_nifti_mrs(tmp_path / 'raw_RAWDATA_19.nii.gz')

    assert img.shape == (1, 1, 1, 2048, 4, 32)
    assert np.iscomplexobj(img.dataobj)
    assert np.isclose(1 / img.header['pixdim'][4], 3012.05)

    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert np.isclose(hdr_ext['SpectrometerFrequency'][0], 300.32)
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == str(data_path.absolute() / '19' / 'rawdata.job0')
    assert 'method' in hdr_ext
    assert 'acqp' in hdr_ext
